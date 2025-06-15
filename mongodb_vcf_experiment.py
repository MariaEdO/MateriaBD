#!/usr/bin/env python3
"""
MongoDB VCF Performance Experiment
Standalone script to load VCF data into MongoDB and measure query performance
"""

import pymongo
from pymongo import MongoClient, ASCENDING, DESCENDING
import time
import statistics
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
from datetime import datetime
from bson import ObjectId

class VCFParser:
    """Custom VCF parser to handle VCF files without external dependencies"""
    
    def __init__(self, vcf_file_path):
        self.vcf_file_path = vcf_file_path
        self.samples = []
        self.header_lines = []
        
    def parse_header(self):
        """Parse VCF header to extract sample names"""
        with open(self.vcf_file_path, 'r') as f:
            for line in f:
                if line.startswith('##'):
                    self.header_lines.append(line.strip())
                elif line.startswith('#CHROM'):
                    # This is the column header line
                    columns = line.strip().split('\t')
                    # Samples start after the first 9 standard columns
                    self.samples = columns[9:] if len(columns) > 9 else []
                    break
    
    def parse_variants(self):
        """Generator that yields parsed variant records"""
        self.parse_header()
        
        with open(self.vcf_file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 8:
                    print(f"Warning: Skipping malformed line {line_num}")
                    continue
                
                try:
                    # Parse basic variant info
                    variant = {
                        'chrom': fields[0],
                        'pos': int(fields[1]) if fields[1] != '.' else None,
                        'id': fields[2] if fields[2] != '.' else None,
                        'ref': fields[3],
                        'alt': fields[4].split(',') if fields[4] != '.' else [],
                        'qual': float(fields[5]) if fields[5] not in ['.', ''] else None,
                        'filter': fields[6].split(';') if fields[6] != '.' else [],
                        'info': self.parse_info(fields[7]) if len(fields) > 7 else {},
                        'format': fields[8].split(':') if len(fields) > 8 else [],
                        'samples': {}
                    }
                    
                    # Parse sample data
                    if len(fields) > 9 and variant['format']:
                        for i, sample_name in enumerate(self.samples):
                            if i + 9 < len(fields):
                                sample_data = fields[i + 9].split(':')
                                variant['samples'][sample_name] = {}
                                for j, fmt_key in enumerate(variant['format']):
                                    if j < len(sample_data) and sample_data[j] != '.':
                                        variant['samples'][sample_name][fmt_key] = sample_data[j]
                    
                    yield variant
                    
                except Exception as e:
                    print(f"Warning: Error parsing line {line_num}: {e}")
                    continue
    
    def parse_info(self, info_string):
        """Parse INFO field into dictionary"""
        info_dict = {}
        if info_string == '.':
            return info_dict
        
        for item in info_string.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                # Try to convert to appropriate type
                if ',' in value:
                    # Multiple values
                    values = value.split(',')
                    try:
                        info_dict[key] = [float(v) for v in values]
                    except ValueError:
                        info_dict[key] = values
                else:
                    # Single value
                    try:
                        info_dict[key] = float(value)
                    except ValueError:
                        info_dict[key] = value
            else:
                # Flag (no value)
                info_dict[item] = True
        
        return info_dict

class MongoDBVCFExperiment:
    """Main experiment class for MongoDB VCF performance testing"""
    
    def __init__(self, vcf_file_path, db_params=None):
        self.vcf_file_path = vcf_file_path
        self.db_params = db_params or {
            'host': 'localhost',
            'port': 27017,
            'database': 'vcf_db'
        }
        self.results = {}
        self.client = None
        self.db = None
        
    def test_connection(self):
        """Test database connection"""
        try:
            self.client = MongoClient(
                host=self.db_params['host'],
                port=self.db_params['port'],
                serverSelectionTimeoutMS=5000
            )
            # Test connection
            self.client.server_info()
            self.db = self.client[self.db_params['database']]
            
            print(f"✅ MongoDB connected successfully!")
            print(f"   Host: {self.db_params['host']}:{self.db_params['port']}")
            print(f"   Database: {self.db_params['database']}")
            return True
        except Exception as e:
            print(f"❌ MongoDB connection failed: {e}")
            print("Make sure MongoDB is running:")
            print("docker run --name mongo-vcf -p 27017:27017 -d mongo:7.0")
            return False
    
    def create_schema(self):
        """Create MongoDB collections and indexes"""
        print("Creating MongoDB collections and indexes...")
        
        # Drop existing collections
        self.db.variants.drop()
        self.db.samples.drop()
        self.db.genotypes.drop()
        
        # Create collections (MongoDB creates them automatically on first insert)
        variants_collection = self.db.variants
        samples_collection = self.db.samples
        genotypes_collection = self.db.genotypes
        
        # We'll create indexes after data loading for better performance
        print("✅ Collections prepared!")
    
    def create_indexes(self):
        """Create indexes after data loading"""
        print("Creating indexes...")
        
        # Variant indexes
        self.db.variants.create_index([("chrom", ASCENDING), ("pos", ASCENDING)])
        self.db.variants.create_index([("qual", ASCENDING)])
        self.db.variants.create_index([("chrom", ASCENDING)])
        self.db.variants.create_index([("info.DP", ASCENDING)])
        
        # Genotype indexes
        self.db.genotypes.create_index([("variant_id", ASCENDING), ("sample_name", ASCENDING)])
        self.db.genotypes.create_index([("sample_name", ASCENDING)])
        self.db.genotypes.create_index([("genotype", ASCENDING)])
        
        # Sample indexes
        self.db.samples.create_index([("sample_name", ASCENDING)], unique=True)
        
        print("✅ Indexes created!")
    
    def load_data(self):
        """Load VCF data into MongoDB"""
        print(f"Loading VCF data from {self.vcf_file_path}...")
        
        if not os.path.exists(self.vcf_file_path):
            raise FileNotFoundError(f"VCF file not found: {self.vcf_file_path}")
        
        # Parse VCF file
        vcf_parser = VCFParser(self.vcf_file_path)
        vcf_parser.parse_header()
        
        print(f"Found {len(vcf_parser.samples)} samples in VCF file")
        
        # Insert samples
        sample_docs = []
        for sample in vcf_parser.samples:
            sample_docs.append({"sample_name": sample})
        
        if sample_docs:
            try:
                self.db.samples.insert_many(sample_docs, ordered=False)
            except Exception as e:
                # Handle duplicate key errors (samples already exist)
                pass
        
        start_time = time.time()
        variant_count = 0
        genotype_count = 0
        
        # Batch inserts for better performance
        variant_batch = []
        genotype_batch = []
        batch_size = 1000
        
        for variant in vcf_parser.parse_variants():
            if variant['pos'] is None:
                continue
            
            # Prepare variant document
            variant_doc = {
                'chrom': variant['chrom'],
                'pos': variant['pos'],
                'variant_id': variant['id'],
                'ref': variant['ref'],
                'alt': variant['alt'],
                'qual': variant['qual'],
                'filter': variant['filter'],
                'info': variant['info'],
                'created_at': datetime.now()
            }
            
            variant_batch.append(variant_doc)
            
            # Prepare genotype documents
            for sample_name, sample_data in variant['samples'].items():
                genotype_doc = {
                    'variant_chrom': variant['chrom'],
                    'variant_pos': variant['pos'],
                    'sample_name': sample_name,
                    'genotype': sample_data.get('GT'),
                    'genotype_quality': int(sample_data.get('GQ', 0)) if sample_data.get('GQ', '').replace('.', '').isdigit() else None,
                    'depth': int(sample_data.get('DP', 0)) if sample_data.get('DP', '').replace('.', '').isdigit() else None,
                    'other_format': sample_data
                }
                genotype_batch.append(genotype_doc)
            
            variant_count += 1
            
            # Insert batches when they reach batch_size
            if len(variant_batch) >= batch_size:
                # Insert variants
                result = self.db.variants.insert_many(variant_batch)
                
                # Add variant ObjectIds to genotype documents
                variant_ids = result.inserted_ids
                genotype_batch_with_ids = []
                
                for i, genotype_doc in enumerate(genotype_batch):
                    # Find the corresponding variant ID
                    variant_idx = i // len(vcf_parser.samples) if vcf_parser.samples else 0
                    if variant_idx < len(variant_ids):
                        genotype_doc['variant_id'] = variant_ids[variant_idx]
                    genotype_batch_with_ids.append(genotype_doc)
                
                # Insert genotypes
                if genotype_batch_with_ids:
                    self.db.genotypes.insert_many(genotype_batch_with_ids)
                    genotype_count += len(genotype_batch_with_ids)
                
                variant_batch = []
                genotype_batch = []
                
                if variant_count % 1000 == 0:
                    print(f"Loaded {variant_count} variants, {genotype_count} genotypes...")
        
        # Insert remaining batches
        if variant_batch:
            result = self.db.variants.insert_many(variant_batch)
            
            if genotype_batch:
                variant_ids = result.inserted_ids
                genotype_batch_with_ids = []
                
                for i, genotype_doc in enumerate(genotype_batch):
                    variant_idx = i // len(vcf_parser.samples) if vcf_parser.samples else 0
                    if variant_idx < len(variant_ids):
                        genotype_doc['variant_id'] = variant_ids[variant_idx]
                    genotype_batch_with_ids.append(genotype_doc)
                
                self.db.genotypes.insert_many(genotype_batch_with_ids)
                genotype_count += len(genotype_batch_with_ids)
        
        load_time = time.time() - start_time
        
        # Create indexes after data loading
        self.create_indexes()
        
        self.results['load_time'] = load_time
        self.results['variant_count'] = variant_count
        self.results['genotype_count'] = genotype_count
        self.results['sample_count'] = len(vcf_parser.samples)
        
        print(f"✅ Data loading completed!")
        print(f"   Loaded {variant_count} variants")
        print(f"   Loaded {genotype_count} genotypes")
        print(f"   Processing {len(vcf_parser.samples)} samples")
        print(f"   Load time: {load_time:.2f} seconds")
        
        return load_time, variant_count
    
    def run_performance_tests(self):
        """Run various query performance tests"""
        print("Running performance tests...")
        
        # Get some sample data for parameterized queries
        sample_variant = self.db.variants.find_one()
        sample_chrom = sample_variant['chrom'] if sample_variant else '1'
        
        sample_sample = self.db.samples.find_one()
        sample_name = sample_sample['sample_name'] if sample_sample else None
        
        # Get position range for the sample chromosome
        pipeline = [
            {"$match": {"chrom": sample_chrom}},
            {"$group": {
                "_id": None,
                "min_pos": {"$min": "$pos"},
                "max_pos": {"$max": "$pos"}
            }}
        ]
        pos_stats = list(self.db.variants.aggregate(pipeline))
        
        if pos_stats:
            min_pos = pos_stats[0]['min_pos']
            max_pos = pos_stats[0]['max_pos']
        else:
            min_pos = 1000000
            max_pos = 2000000
        
        queries = {
            "total_variants": {
                "operation": lambda: self.db.variants.count_documents({}),
                "description": "Count total variants"
            },
            "range_query": {
                "operation": lambda: self.db.variants.count_documents({
                    "chrom": sample_chrom,
                    "pos": {"$gte": min_pos, "$lte": min_pos + 100000}
                }),
                "description": f"Count variants in range {sample_chrom}:{min_pos}-{min_pos + 100000}"
            },
            "quality_filter": {
                "operation": lambda: self.db.variants.count_documents({"qual": {"$gt": 30}}),
                "description": "Count high-quality variants (QUAL > 30)"
            },
            "chromosome_summary": {
                "operation": lambda: list(self.db.variants.aggregate([
                    {"$group": {
                        "_id": "$chrom",
                        "variant_count": {"$sum": 1},
                        "avg_qual": {"$avg": "$qual"}
                    }},
                    {"$sort": {"_id": 1}}
                ])),
                "description": "Variants per chromosome with average quality"
            },
            "info_field_query": {
                "operation": lambda: self.db.variants.count_documents({"info.DP": {"$exists": True}}),
                "description": "Count variants with depth information"
            },
            "complex_aggregation": {
                "operation": lambda: list(self.db.variants.aggregate([
                    {"$match": {"chrom": sample_chrom}},
                    {"$lookup": {
                        "from": "genotypes",
                        "localField": "_id",
                        "foreignField": "variant_id",
                        "as": "genotypes"
                    }},
                    {"$project": {
                        "chrom": 1,
                        "pos": 1,
                        "sample_count": {"$size": "$genotypes"},
                        "avg_depth": {"$avg": "$genotypes.depth"}
                    }},
                    {"$sort": {"pos": 1}},
                    {"$limit": 100}
                ])),
                "description": f"Complex aggregation for chromosome {sample_chrom}"
            }
        }
        
        # Add sample-specific query if we have samples
        if sample_name:
            queries["sample_genotypes"] = {
                "operation": lambda: list(self.db.genotypes.aggregate([
                    {"$match": {"sample_name": sample_name}},
                    {"$lookup": {
                        "from": "variants",
                        "localField": "variant_id",
                        "foreignField": "_id",
                        "as": "variant"
                    }},
                    {"$unwind": "$variant"},
                    {"$project": {
                        "chrom": "$variant.chrom",
                        "pos": "$variant.pos",
                        "genotype": 1,
                        "genotype_quality": 1
                    }},
                    {"$limit": 100}
                ])),
                "description": f"Genotypes for sample {sample_name}"
            }
        
        query_results = {}
        
        for query_name, query_info in queries.items():
            times = []
            result_count = 0
            
            print(f"  Testing: {query_info['description']}")
            
            # Run each query 3 times for better average
            for run in range(3):
                start_time = time.time()
                
                try:
                    result = query_info['operation']()
                    query_time = time.time() - start_time
                    times.append(query_time)
                    
                    if run == 0:  # Store result count from first run
                        if isinstance(result, list):
                            result_count = len(result)
                        else:
                            result_count = result  # For count queries
                        
                except Exception as e:
                    print(f"    Error in query {query_name}: {e}")
                    times.append(float('inf'))
            
            # Filter out failed queries
            valid_times = [t for t in times if t != float('inf')]
            
            if valid_times:
                query_results[query_name] = {
                    'description': query_info['description'],
                    'avg_time': statistics.mean(valid_times),
                    'min_time': min(valid_times),
                    'max_time': max(valid_times),
                    'result_count': result_count,
                    'times': valid_times
                }
                print(f"    Average time: {statistics.mean(valid_times)*1000:.2f} ms")
            else:
                print(f"    Query failed on all attempts")
        
        self.results['queries'] = query_results
        return query_results
    
    def get_database_statistics(self):
        """Get database size and statistics"""
        print("Collecting database statistics...")
        
        stats = {}
        
        # Database statistics
        db_stats = self.db.command("dbStats")
        stats['database_size'] = f"{db_stats['dataSize'] / (1024*1024):.2f} MB"
        stats['index_size'] = f"{db_stats['indexSize'] / (1024*1024):.2f} MB"
        stats['total_size'] = f"{db_stats['storageSize'] / (1024*1024):.2f} MB"
        
        # Collection statistics
        collection_stats = []
        for collection_name in ['variants', 'samples', 'genotypes']:
            try:
                coll_stats = self.db.command("collStats", collection_name)
                collection_stats.append({
                    'name': collection_name,
                    'count': coll_stats['count'],
                    'size': f"{coll_stats['size'] / (1024*1024):.2f} MB",
                    'avgObjSize': f"{coll_stats.get('avgObjSize', 0):.2f} bytes"
                })
            except Exception as e:
                print(f"Could not get stats for collection {collection_name}: {e}")
        
        stats['collection_stats'] = collection_stats
        
        # Index statistics
        index_stats = []
        for collection_name in ['variants', 'samples', 'genotypes']:
            try:
                indexes = self.db[collection_name].list_indexes()
                for index in indexes:
                    index_stats.append({
                        'collection': collection_name,
                        'name': index.get('name', 'unknown'),
                        'keys': str(index.get('key', {}))
                    })
            except Exception as e:
                print(f"Could not get index stats for {collection_name}: {e}")
        
        stats['index_stats'] = index_stats
        
        self.results['database_stats'] = stats
        print(f"✅ Database size: {stats['database_size']}")
        
        return stats
    
    def generate_report(self, output_file="mongodb_vcf_report.txt"):
        """Generate a detailed report of the experiment"""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        report = f"""
MongoDB VCF Performance Experiment Report
Generated: {timestamp}
VCF File: {self.vcf_file_path}

=== DATA LOADING RESULTS ===
Load Time: {self.results.get('load_time', 0):.2f} seconds
Variants Loaded: {self.results.get('variant_count', 0):,}
Genotypes Loaded: {self.results.get('genotype_count', 0):,}
Samples: {self.results.get('sample_count', 0)}
Loading Rate: {self.results.get('variant_count', 0) / max(self.results.get('load_time', 1), 1):.0f} variants/second

=== DATABASE STATISTICS ===
"""
        
        if 'database_stats' in self.results:
            stats = self.results['database_stats']
            report += f"Database Size: {stats.get('database_size', 'Unknown')}\n"
            report += f"Index Size: {stats.get('index_size', 'Unknown')}\n"
            report += f"Total Size: {stats.get('total_size', 'Unknown')}\n\n"
            
            report += "Collection Statistics:\n"
            for coll_stat in stats.get('collection_stats', []):
                report += f"  {coll_stat['name']}: {coll_stat['count']:,} documents, {coll_stat['size']}\n"
            
            report += "\nIndexes:\n"
            for index_stat in stats.get('index_stats', []):
                report += f"  {index_stat['collection']}.{index_stat['name']}: {index_stat['keys']}\n"
        
        report += "\n=== QUERY PERFORMANCE ===\n"
        
        if 'queries' in self.results:
            for query_name, query_data in self.results['queries'].items():
                report += f"""
{query_data['description']}:
  Average Time: {query_data['avg_time']*1000:.2f} ms
  Min Time: {query_data['min_time']*1000:.2f} ms
  Max Time: {query_data['max_time']*1000:.2f} ms
  Results Returned: {query_data['result_count']:,}
"""
        
        # Write to file
        with open(output_file, 'w') as f:
            f.write(report)
        
        print(f"✅ Report saved to {output_file}")
        print(report)
        
        return report
    
    def create_visualizations(self):
        """Create performance visualization charts"""
        if 'queries' not in self.results:
            print("No query results to visualize")
            return
        
        queries = self.results['queries']
        
        # Prepare data for plotting
        query_names = list(queries.keys())
        avg_times = [queries[name]['avg_time'] * 1000 for name in query_names]  # Convert to ms
        
        # Create figure with subplots
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
        
        # Query performance bar chart
        bars = ax1.bar(range(len(query_names)), avg_times, alpha=0.7, color='forestgreen')
        ax1.set_xlabel('Query Type')
        ax1.set_ylabel('Average Time (ms)')
        ax1.set_title('MongoDB Query Performance')
        ax1.set_xticks(range(len(query_names)))
        ax1.set_xticklabels([name.replace('_', '\n') for name in query_names], rotation=45, ha='right')
        ax1.set_yscale('log')  # Log scale for better visibility
        
        # Add value labels on bars
        for bar, time_val in zip(bars, avg_times):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{time_val:.1f}ms',
                    ha='center', va='bottom', fontsize=8)
        
        # Data loading and summary statistics
        if all(key in self.results for key in ['load_time', 'variant_count', 'sample_count']):
            stats_labels = ['Load Time (s)', 'Variants (k)', 'Samples', 'Genotypes (k)']
            stats_values = [
                self.results['load_time'],
                self.results['variant_count'] / 1000,
                self.results['sample_count'],
                self.results['genotype_count'] / 1000
            ]
            
            bars2 = ax2.bar(stats_labels, stats_values, alpha=0.7, color=['orange', 'green', 'red', 'purple'])
            ax2.set_ylabel('Value')
            ax2.set_title('Data Loading Statistics')
            ax2.set_yscale('log')
            
            # Add value labels
            for bar, val in zip(bars2, stats_values):
                height = bar.get_height()
                ax2.text(bar.get_x() + bar.get_width()/2., height,
                        f'{val:.1f}',
                        ha='center', va='bottom')
        
        plt.tight_layout()
        plt.savefig('mongodb_vcf_performance.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        print("✅ Visualization saved as 'mongodb_vcf_performance.png'")
    
    def cleanup(self):
        """Close database connection"""
        if self.client:
            self.client.close()
    
    def run_experiment(self):
        """Run the complete experiment"""
        print("=== MongoDB VCF Performance Experiment ===\n")
        
        # Test connection
        if not self.test_connection():
            return False
        
        try:
            # Create schema
            self.create_schema()
            
            # Load data
            self.load_data()
            
            # Run performance tests
            self.run_performance_tests()
            
            # Get database statistics
            self.get_database_statistics()
            
            # Generate report
            self.generate_report()
            
            # Create visualizations
            self.create_visualizations()
            
            print("\n✅ Experiment completed successfully!")
            return True
            
        except Exception as e:
            print(f"❌ Experiment failed: {e}")
            import traceback
            traceback.print_exc()
            return False
        finally:
            self.cleanup()

def main():
    """Main function"""
    if len(sys.argv) != 2:
        print("Usage: python mongodb_vcf_experiment.py <vcf_file_path>")
        print("Example: python mongodb_vcf_experiment.py test.vcf")
        sys.exit(1)
    
    vcf_file = sys.argv[1]
    
    if not os.path.exists(vcf_file):
        print(f"Error: VCF file '{vcf_file}' not found")
        sys.exit(1)
    
    # Create and run experiment
    experiment = MongoDBVCFExperiment(vcf_file)
    success = experiment.run_experiment()
    
    if success:
        print("\nExperiment Results Summary:")
        print(f"- VCF file processed: {vcf_file}")
        print(f"- Variants loaded: {experiment.results.get('variant_count', 0):,}")
        print(f"- Load time: {experiment.results.get('load_time', 0):.2f} seconds")
        print(f"- Report saved: mongodb_vcf_report.txt")
        print(f"- Chart saved: mongodb_vcf_performance.png")
    else:
        print("Experiment failed. Check the error messages above.")
        sys.exit(1)

if __name__ == "__main__":
    main()
