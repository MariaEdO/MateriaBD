#!/usr/bin/env python3
"""
Enhanced MongoDB VCF Performance Experiment
Standalone script to load VCF data into MongoDB and measure query performance
with comprehensive system resource monitoring
"""

import pymongo
from pymongo import MongoClient, InsertOne, UpdateOne
import time
import statistics
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import threading
import psutil
import json
from datetime import datetime
from collections import defaultdict
import math
from bson import ObjectId

def clean_for_json(obj):
    """
    Recursively clean an object to make it JSON-serializable by handling NaN, Infinity, etc.
    """
    if isinstance(obj, dict):
        return {k: clean_for_json(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [clean_for_json(item) for item in obj]
    elif isinstance(obj, float):
        if math.isnan(obj):
            return None
        elif math.isinf(obj):
            return None
        else:
            return obj
    elif isinstance(obj, ObjectId):
        return str(obj)
    else:
        return obj

class SystemMonitor:
    """Monitor system resources during the experiment"""
    
    def __init__(self, interval=1.0):
        self.interval = interval
        self.monitoring = False
        self.metrics = {
            'timestamps': [],
            'cpu_percent': [],
            'memory_percent': [],
            'memory_used_gb': [],
            'disk_read_mb': [],
            'disk_write_mb': [],
            'network_sent_mb': [],
            'network_recv_mb': []
        }
        self.process = psutil.Process(os.getpid())
        
        # Initialize disk and network counters with error handling
        try:
            self.initial_disk_io = psutil.disk_io_counters()
            self.initial_network_io = psutil.net_io_counters()
        except Exception as e:
            print(f"Warning: Could not initialize I/O counters: {e}")
            self.initial_disk_io = None
            self.initial_network_io = None
        
    def start_monitoring(self):
        """Start monitoring system resources"""
        self.monitoring = True
        self.monitor_thread = threading.Thread(target=self._monitor_loop)
        self.monitor_thread.daemon = True
        self.monitor_thread.start()
        print("🔍 System monitoring started")
        
    def stop_monitoring(self):
        """Stop monitoring system resources"""
        self.monitoring = False
        if hasattr(self, 'monitor_thread'):
            self.monitor_thread.join(timeout=2)
        print("🔍 System monitoring stopped")
        
    def _monitor_loop(self):
        """Main monitoring loop"""
        while self.monitoring:
            try:
                timestamp = time.time()
                
                # CPU and memory usage
                cpu_percent = psutil.cpu_percent(interval=None)
                memory = psutil.virtual_memory()
                
                # Process-specific memory
                process_memory = self.process.memory_info()
                
                # Disk I/O with error handling
                disk_read_mb = 0
                disk_write_mb = 0
                if self.initial_disk_io:
                    try:
                        current_disk_io = psutil.disk_io_counters()
                        if current_disk_io:
                            disk_read_mb = (current_disk_io.read_bytes - self.initial_disk_io.read_bytes) / (1024 * 1024)
                            disk_write_mb = (current_disk_io.write_bytes - self.initial_disk_io.write_bytes) / (1024 * 1024)
                    except Exception:
                        pass
                
                # Network I/O with error handling
                network_sent_mb = 0
                network_recv_mb = 0
                if self.initial_network_io:
                    try:
                        current_network_io = psutil.net_io_counters()
                        if current_network_io:
                            network_sent_mb = (current_network_io.bytes_sent - self.initial_network_io.bytes_sent) / (1024 * 1024)
                            network_recv_mb = (current_network_io.bytes_recv - self.initial_network_io.bytes_recv) / (1024 * 1024)
                    except Exception:
                        pass
                
                # Store metrics
                self.metrics['timestamps'].append(timestamp)
                self.metrics['cpu_percent'].append(cpu_percent)
                self.metrics['memory_percent'].append(memory.percent)
                self.metrics['memory_used_gb'].append(process_memory.rss / (1024**3))
                self.metrics['disk_read_mb'].append(disk_read_mb)
                self.metrics['disk_write_mb'].append(disk_write_mb)
                self.metrics['network_sent_mb'].append(network_sent_mb)
                self.metrics['network_recv_mb'].append(network_recv_mb)
                
                time.sleep(self.interval)
                
            except Exception as e:
                print(f"Warning: Error in system monitoring: {e}")
                time.sleep(self.interval)
    
    def get_summary_stats(self):
        """Get summary statistics for all monitored metrics"""
        if not self.metrics['timestamps']:
            return {}
            
        stats = {}
        for metric_name, values in self.metrics.items():
            if metric_name == 'timestamps':
                continue
            if values:
                stats[metric_name] = {
                    'avg': statistics.mean(values),
                    'max': max(values),
                    'min': min(values),
                    'final': values[-1] if values else 0
                }
        
        # Add duration
        if len(self.metrics['timestamps']) > 1:
            stats['duration_seconds'] = self.metrics['timestamps'][-1] - self.metrics['timestamps'][0]
            
        return stats

class VCFParser:
    """Custom VCF parser to handle VCF files without external dependencies"""
    
    def __init__(self, vcf_file_path):
        self.vcf_file_path = vcf_file_path
        self.samples = []
        self.header_lines = []
        self.file_stats = {}
        
    def analyze_file(self):
        """Analyze VCF file and collect basic statistics"""
        print("🔍 Analyzing VCF file...")
        
        file_size = os.path.getsize(self.vcf_file_path)
        self.file_stats['file_size_mb'] = file_size / (1024 * 1024)
        
        line_count = 0
        header_lines = 0
        variant_lines = 0
        
        with open(self.vcf_file_path, 'r') as f:
            for line in f:
                line_count += 1
                if line.startswith('#'):
                    header_lines += 1
                else:
                    variant_lines += 1
                    
        self.file_stats.update({
            'total_lines': line_count,
            'header_lines': header_lines,
            'variant_lines': variant_lines
        })
        
        print(f"   File size: {self.file_stats['file_size_mb']:.2f} MB")
        print(f"   Total lines: {line_count:,}")
        print(f"   Variant lines: {variant_lines:,}")
        
        return self.file_stats
        
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
                        'variant_id': fields[2] if fields[2] != '.' else None,
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
                                        # Convert numeric values
                                        value = sample_data[j]
                                        if fmt_key in ['GQ', 'DP'] and value.replace('.', '').isdigit():
                                            try:
                                                variant['samples'][sample_name][fmt_key] = int(float(value))
                                            except ValueError:
                                                variant['samples'][sample_name][fmt_key] = value
                                        else:
                                            variant['samples'][sample_name][fmt_key] = value
                    
                    yield variant
                    
                except Exception as e:
                    print(f"Warning: Error parsing line {line_num}: {e}")
                    continue

    def parse_info(self, info_string):
        """Parse INFO field into dictionary, handling NaN and other special values"""
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
                    converted_values = []
                    for v in values:
                        try:
                            # Try float first
                            float_val = float(v)
                            converted_values.append(float_val)
                        except ValueError:
                            # Keep as string if not numeric
                            converted_values.append(v)
                    info_dict[key] = converted_values
                else:
                    # Single value
                    try:
                        # Try float first
                        float_val = float(value)
                        info_dict[key] = float_val
                    except ValueError:
                        # Keep as string if not numeric
                        info_dict[key] = value
            else:
                # Flag (no value)
                info_dict[item] = True
    
        # Clean the dictionary to handle NaN, Infinity, etc.
        return clean_for_json(info_dict)

class MongoDBVCFExperiment:
    """Main experiment class for MongoDB VCF performance testing"""
    
    def __init__(self, vcf_file_path, mongo_params=None):
        self.vcf_file_path = vcf_file_path
        self.mongo_params = mongo_params or {
            'host': 'localhost',
            'port': 27017,
            'database': 'vcf_db'
        }
        self.results = {}
        self.system_monitor = SystemMonitor()
        self.client = None
        self.db = None
        
    def test_connection(self):
        """Test MongoDB connection"""
        try:
            # Create connection with timeout
            self.client = MongoClient(
                host=self.mongo_params['host'],
                port=self.mongo_params['port'],
                serverSelectionTimeoutMS=5000
            )
            
            # Test connection
            self.client.admin.command('ping')
            
            # Get database
            self.db = self.client[self.mongo_params['database']]
            
            # Get server info
            server_info = self.client.server_info()
            print(f"✅ MongoDB connected successfully!")
            print(f"   Version: {server_info['version']}")
            print(f"   Database: {self.mongo_params['database']}")
            
            return True
            
        except Exception as e:
            print(f"❌ MongoDB connection failed: {e}")
            print("Make sure MongoDB is running:")
            print("docker run --name mongo-vcf -p 27017:27017 -d mongo:latest")
            return False
    
    def create_schema(self):
        """Create MongoDB collections and indexes"""
        print("Creating MongoDB schema...")
        
        start_time = time.time()
        
        try:
            # Drop existing collections
            self.db.variants.drop()
            self.db.samples.drop()
            
            # Create collections (MongoDB creates them automatically on first insert)
            # But we can create indexes
            
            # Indexes for variants collection
            variant_indexes = [
                [("chrom", 1), ("pos", 1)],  # Compound index for range queries
                [("chrom", 1)],              # Chromosome queries
                [("pos", 1)],                # Position queries
                [("qual", 1)],               # Quality filtering
                [("variant_id", 1)],         # ID lookups
                [("info.DP", 1)],            # INFO field queries
                [("samples.*.GT", 1)],       # Genotype queries (wildcard)
            ]
            
            for index_spec in variant_indexes:
                try:
                    self.db.variants.create_index(index_spec)
                except Exception as e:
                    print(f"Warning: Could not create index {index_spec}: {e}")
            
            # Indexes for samples collection
            self.db.samples.create_index([("sample_name", 1)], unique=True)
            
            # Text indexes for flexible searching
            try:
                self.db.variants.create_index([("$**", "text")])
            except Exception as e:
                print(f"Warning: Could not create text index: {e}")
            
            print("✅ MongoDB collections and indexes created")
            
        except Exception as e:
            print(f"Warning: Error in schema creation: {e}")
        
        schema_time = time.time() - start_time
        self.results['schema_creation_time'] = schema_time
        
        print(f"✅ Schema created successfully! ({schema_time:.2f}s)")
    
    def load_data(self):
        """Load VCF data into MongoDB with resource monitoring"""
        print(f"Loading VCF data from {self.vcf_file_path}...")
        
        if not os.path.exists(self.vcf_file_path):
            raise FileNotFoundError(f"VCF file not found: {self.vcf_file_path}")
        
        # Analyze file first
        vcf_parser = VCFParser(self.vcf_file_path)
        file_stats = vcf_parser.analyze_file()
        self.results['file_stats'] = file_stats
        
        # Start system monitoring
        self.system_monitor.start_monitoring()
        
        try:
            # Parse VCF file
            vcf_parser.parse_header()
            
            print(f"Found {len(vcf_parser.samples)} samples in VCF file")
            
            # Insert samples first
            sample_insert_start = time.time()
            sample_docs = []
            for i, sample in enumerate(vcf_parser.samples):
                sample_docs.append({
                    '_id': i + 1,  # Use integer IDs for easier referencing
                    'sample_name': sample,
                    'created_at': datetime.utcnow()
                })
            
            if sample_docs:
                self.db.samples.insert_many(sample_docs, ordered=False)
            
            sample_insert_time = time.time() - sample_insert_start
            
            # Create sample name to ID mapping
            sample_id_map = {doc['sample_name']: doc['_id'] for doc in sample_docs}
            
            # Track loading progress
            start_time = time.time()
            variant_count = 0
            batch_size = 1000
            progress_interval = 5000
            
            # Statistics tracking
            chrom_stats = defaultdict(int)
            qual_distribution = []
            genotype_stats = defaultdict(int)
            
            variant_batch = []
            
            # Process variants
            for variant in vcf_parser.parse_variants():
                if variant['pos'] is None:
                    continue
                
                # Collect statistics
                chrom_stats[variant['chrom']] += 1
                if variant['qual'] is not None:
                    qual_distribution.append(variant['qual'])
                
                # Convert samples data for MongoDB
                mongo_samples = {}
                genotype_count = 0
                
                for sample_name, sample_data in variant['samples'].items():
                    if sample_name in sample_id_map:
                        # Collect genotype statistics
                        gt = sample_data.get('GT', 'Unknown')
                        genotype_stats[gt] += 1
                        genotype_count += 1
                        
                        # Store sample data with sample ID reference
                        mongo_samples[str(sample_id_map[sample_name])] = sample_data
                
                # Create MongoDB document
                variant_doc = {
                    'chrom': variant['chrom'],
                    'pos': variant['pos'],
                    'variant_id': variant['variant_id'],
                    'ref': variant['ref'],
                    'alt': variant['alt'],
                    'qual': variant['qual'],
                    'filter': variant['filter'],
                    'info': variant['info'],
                    'format': variant['format'],
                    'samples': mongo_samples,
                    'sample_count': genotype_count,
                    'created_at': datetime.utcnow()
                }
                
                variant_batch.append(variant_doc)
                variant_count += 1
                
                # Batch insert when batch is full
                if len(variant_batch) >= batch_size:
                    try:
                        result = self.db.variants.insert_many(variant_batch, ordered=False)
                        if variant_count % progress_interval == 0:
                            print(f"   Inserted batch of {len(result.inserted_ids)} variants")
                    except Exception as e:
                        print(f"   Warning: Batch insert error: {e}")
                    
                    variant_batch = []
                    
                    if variant_count % progress_interval == 0:
                        current_time = time.time()
                        elapsed = current_time - start_time
                        rate = variant_count / elapsed
                        print(f"   Progress: {variant_count:,} variants loaded ({rate:.0f} variants/sec)")
            
            # Insert remaining variants
            if variant_batch:
                try:
                    result = self.db.variants.insert_many(variant_batch, ordered=False)
                    print(f"   Inserted final batch of {len(result.inserted_ids)} variants")
                except Exception as e:
                    print(f"   Warning: Final batch insert error: {e}")
            
            # Calculate total genotypes
            total_genotypes = sum(genotype_stats.values())
            
        except Exception as e:
            print(f"Error during data loading: {e}")
            raise
        finally:
            # Stop monitoring and collect stats
            self.system_monitor.stop_monitoring()
        
        load_time = time.time() - start_time
        
        # Store comprehensive results
        self.results.update({
            'load_time': load_time,
            'sample_insert_time': sample_insert_time,
            'variant_count': variant_count,
            'genotype_count': total_genotypes,
            'sample_count': len(vcf_parser.samples),
            'chromosome_distribution': dict(chrom_stats),
            'quality_stats': {
                'count': len(qual_distribution),
                'mean': statistics.mean(qual_distribution) if qual_distribution else 0,
                'median': statistics.median(qual_distribution) if qual_distribution else 0,
                'min': min(qual_distribution) if qual_distribution else 0,
                'max': max(qual_distribution) if qual_distribution else 0
            },
            'genotype_distribution': dict(genotype_stats),
            'loading_rate_variants_per_sec': variant_count / load_time if load_time > 0 else 0,
            'system_stats_loading': self.system_monitor.get_summary_stats()
        })
        
        print(f"✅ Data loading completed!")
        print(f"   Loaded {variant_count:,} variants")
        print(f"   Loaded {total_genotypes:,} genotypes")
        print(f"   Processing {len(vcf_parser.samples)} samples")
        print(f"   Total load time: {load_time:.2f} seconds")
        print(f"   Loading rate: {variant_count / load_time:.0f} variants/second")
        
        return load_time, variant_count
    
    def run_performance_tests(self):
        """Run various query performance tests with resource monitoring"""
        print("Running performance tests...")
        
        # Start fresh monitoring for queries
        query_monitor = SystemMonitor(interval=0.5)
        query_monitor.start_monitoring()
        
        try:
            # Get database statistics first
            total_variants = self.db.variants.count_documents({})
            total_samples = self.db.samples.count_documents({})
            
            # Calculate total genotypes by summing sample_count from all variants
            pipeline = [{"$group": {"_id": None, "total_genotypes": {"$sum": "$sample_count"}}}]
            genotype_result = list(self.db.variants.aggregate(pipeline))
            total_genotypes = genotype_result[0]['total_genotypes'] if genotype_result else 0
            
            # Get some sample data for parameterized queries
            sample_variant = self.db.variants.find_one({})
            sample_chrom = sample_variant['chrom'] if sample_variant else '1'
            
            sample_doc = self.db.samples.find_one({})
            sample_name = sample_doc['sample_name'] if sample_doc else None
            sample_id = str(sample_doc['_id']) if sample_doc else None
            
            # Get position range for the sample chromosome
            chrom_variants = list(self.db.variants.find(
                {"chrom": sample_chrom}, 
                {"pos": 1}
            ).sort("pos", 1).limit(2))
            
            min_pos = chrom_variants[0]['pos'] if chrom_variants else 1000000
            max_pos = chrom_variants[-1]['pos'] if len(chrom_variants) > 1 else min_pos + 100000
            
            # Define MongoDB queries equivalent to PostgreSQL ones
            queries = {
                "total_variants": {
                    "operation": lambda: self.db.variants.count_documents({}),
                    "description": "Count total variants",
                    "expected_result_type": "count"
                },
                "total_genotypes": {
                    "operation": lambda: list(self.db.variants.aggregate([
                        {"$group": {"_id": None, "total": {"$sum": "$sample_count"}}}
                    ]))[0]['total'],
                    "description": "Count total genotypes",
                    "expected_result_type": "count"
                },
                "range_query": {
                    "operation": lambda: self.db.variants.count_documents({
                        "chrom": sample_chrom,
                        "pos": {"$gte": min_pos, "$lte": min_pos + 100000}
                    }),
                    "description": f"Count variants in range {sample_chrom}:{min_pos}-{min_pos + 100000}",
                    "expected_result_type": "count"
                },
                "quality_filter": {
                    "operation": lambda: self.db.variants.count_documents({"qual": {"$gt": 30}}),
                    "description": "Count high-quality variants (QUAL > 30)",
                    "expected_result_type": "count"
                },
                "quality_stats": {
                    "operation": lambda: list(self.db.variants.aggregate([
                        {"$match": {"qual": {"$ne": None}}},
                        {"$group": {
                            "_id": None,
                            "avg_qual": {"$avg": "$qual"},
                            "min_qual": {"$min": "$qual"},
                            "max_qual": {"$max": "$qual"},
                            "count": {"$sum": 1}
                        }}
                    ])),
                    "description": "Quality statistics",
                    "expected_result_type": "stats"
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
                    "description": "Variants per chromosome with average quality",
                    "expected_result_type": "summary"
                },
                "info_field_query": {
                    "operation": lambda: self.db.variants.count_documents({"info.DP": {"$exists": True}}),
                    "description": "Count variants with depth information",
                    "expected_result_type": "count"
                },
                "genotype_distribution": {
                    "operation": lambda: list(self.db.variants.aggregate([
                        {"$project": {
                            "sample_data": {"$objectToArray": "$samples"}
                        }},
                        {"$unwind": "$sample_data"},
                        {"$project": {
                            "genotype": "$sample_data.v.GT"
                        }},
                        {"$match": {"genotype": {"$ne": None, "$exists": True}}},
                        {"$group": {
                            "_id": "$genotype",
                            "count": {"$sum": 1}
                        }},
                        {"$sort": {"count": -1}},
                        {"$limit": 10}
                    ])),
                    "description": "Top genotype distribution",
                    "expected_result_type": "distribution"
                },
                "complex_join": {
                    "operation": lambda: list(self.db.variants.aggregate([
                        {"$match": {"chrom": sample_chrom}},
                        {"$project": {
                            "chrom": 1,
                            "pos": 1,
                            "sample_count": 1,
                            "avg_depth": {
                                "$avg": {
                                    "$map": {
                                        "input": {"$objectToArray": "$samples"},
                                        "as": "sample",
                                        "in": {"$toDouble": "$$sample.v.DP"}
                                    }
                                }
                            }
                        }},
                        {"$sort": {"pos": 1}},
                        {"$limit": 100}
                    ])),
                    "description": f"Complex aggregation for chromosome {sample_chrom}",
                    "expected_result_type": "complex"
                },
                "depth_statistics": {
                    "operation": lambda: list(self.db.variants.aggregate([
                        {"$project": {
                            "depths": {
                                "$map": {
                                    "input": {"$objectToArray": "$samples"},
                                    "as": "sample",
                                    "in": {"$toDouble": "$$sample.v.DP"}
                                }
                            }
                        }},
                        {"$unwind": "$depths"},
                        {"$match": {"depths": {"$ne": None}}},
                        {"$group": {
                            "_id": None,
                            "avg_depth": {"$avg": "$depths"},
                            "min_depth": {"$min": "$depths"},
                            "max_depth": {"$max": "$depths"},
                            "count": {"$sum": 1}
                        }}
                    ])),
                    "description": "Depth statistics across all genotypes",
                    "expected_result_type": "stats"
                }
            }
            
            # Add sample-specific query if we have samples
            if sample_id:
                queries["sample_genotypes"] = {
                    "operation": lambda: list(self.db.variants.find(
                        {f"samples.{sample_id}": {"$exists": True}},
                        {"chrom": 1, "pos": 1, f"samples.{sample_id}": 1}
                    ).limit(100)),
                    "description": f"Genotypes for sample {sample_name}",
                    "expected_result_type": "sample_data"
                }
            
            query_results = {}
            
            for query_name, query_info in queries.items():
                times = []
                result_count = 0
                result_data = None
                
                print(f"  Testing: {query_info['description']}")
                
                # Run each query 5 times for better statistics
                for run in range(5):
                    start_time = time.time()
                    
                    try:
                        result = query_info['operation']()
                        query_time = time.time() - start_time
                        times.append(query_time)
                        
                        if run == 0:  # Store result from first run
                            if isinstance(result, list):
                                result_count = len(result)
                                result_data = result[:10] if result else []
                            else:
                                result_count = 1
                                result_data = [result]
                            
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
                        'std_time': statistics.stdev(valid_times) if len(valid_times) > 1 else 0,
                        'result_count': result_count,
                        'sample_results': clean_for_json(result_data),
                        'times': valid_times,
                        'query_type': query_info['expected_result_type']
                    }
                    print(f"    Average time: {statistics.mean(valid_times)*1000:.2f} ms")
                    print(f"    Results returned: {result_count:,}")
                else:
                    print(f"    Query failed on all attempts")
                    
        except Exception as e:
            print(f"Error during performance tests: {e}")
            raise
        finally:
            # Stop query monitoring
            query_monitor.stop_monitoring()
        
        self.results['queries'] = query_results
        self.results['system_stats_queries'] = query_monitor.get_summary_stats()
        self.results['database_totals'] = {
            'total_variants': total_variants,
            'total_genotypes': total_genotypes,
            'total_samples': total_samples
        }
        
        return query_results
    
    def get_database_statistics(self):
        """Get comprehensive database size and statistics"""
        print("Collecting database statistics...")
        
        try:
            stats = {}
            
            # Database statistics
            db_stats = self.db.command("dbStats")
            stats['database_size_bytes'] = db_stats.get('dataSize', 0)
            stats['database_size_mb'] = stats['database_size_bytes'] / (1024 * 1024)
            stats['index_size_bytes'] = db_stats.get('indexSize', 0)
            stats['index_size_mb'] = stats['index_size_bytes'] / (1024 * 1024)
            stats['storage_size_bytes'] = db_stats.get('storageSize', 0)
            stats['storage_size_mb'] = stats['storage_size_bytes'] / (1024 * 1024)
            stats['collection_count'] = db_stats.get('collections', 0)
            stats['object_count'] = db_stats.get('objects', 0)
            
            # Collection-specific statistics
            collection_stats = {}
            
            for collection_name in ['variants', 'samples']:
                try:
                    coll_stats = self.db.command("collStats", collection_name)
                    collection_stats[collection_name] = {
                        'count': coll_stats.get('count', 0),
                        'size_bytes': coll_stats.get('size', 0),
                        'size_mb': coll_stats.get('size', 0) / (1024 * 1024),
                        'storage_size_bytes': coll_stats.get('storageSize', 0),
                        'storage_size_mb': coll_stats.get('storageSize', 0) / (1024 * 1024),
                        'avg_obj_size': coll_stats.get('avgObjSize', 0),
                        'index_count': coll_stats.get('nindexes', 0),
                        'total_index_size_bytes': coll_stats.get('totalIndexSize', 0),
                        'total_index_size_mb': coll_stats.get('totalIndexSize', 0) / (1024 * 1024)
                    }
                except Exception as e:
                    print(f"Warning: Could not get stats for collection {collection_name}: {e}")
                    collection_stats[collection_name] = {}
            
            stats['collection_stats'] = collection_stats
            
            # Index information
            try:
                variant_indexes = list(self.db.variants.list_indexes())
                sample_indexes = list(self.db.samples.list_indexes())
                
                stats['index_info'] = {
                    'variants': [clean_for_json(idx) for idx in variant_indexes],
                    'samples': [clean_for_json(idx) for idx in sample_indexes]
                }
            except Exception as e:
                print(f"Warning: Could not get index information: {e}")
                stats['index_info'] = {}
            
            # Server status
            try:
                server_status = self.db.command("serverStatus")
                stats['server_info'] = {
                    'version': server_status.get('version', 'Unknown'),
                    'uptime': server_status.get('uptime', 0),
                    'connections': server_status.get('connections', {}),
                    'memory': server_status.get('mem', {}),
                    'opcounters': server_status.get('opcounters', {})
                }
            except Exception as e:
                print(f"Warning: Could not get server status: {e}")
                stats['server_info'] = {}
            
            # Calculate some derived metrics
            if stats['database_size_bytes'] > 0:
                variant_count = self.results.get('variant_count', 0)
                stats['storage_efficiency'] = {
                    'variants_per_mb': variant_count / max(stats['database_size_mb'], 1),
                    'bytes_per_variant': stats['database_size_bytes'] / max(variant_count, 1)
                }
                
        except Exception as e:
            print(f"Warning: Error collecting database statistics: {e}")
            stats = {}
        
        self.results['database_stats'] = stats
        
        # Print summary
        if stats:
            print(f"✅ Database size: {stats.get('database_size_mb', 0):.2f} MB")
            print(f"   Storage size: {stats.get('storage_size_mb', 0):.2f} MB")
            print(f"   Index size: {stats.get('index_size_mb', 0):.2f} MB")
            print(f"   Total objects: {stats.get('object_count', 0):,}")
        
        return stats
    
    def generate_comprehensive_report(self, output_file=None):
        """Generate a detailed comprehensive report of the experiment"""
        if output_file is None:
            file_prefix = self._get_file_prefix()
            output_file = f"{file_prefix}_mongodb_vcf_comprehensive_report.txt"
            
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        report = f"""
==============================================================================
MongoDB VCF Performance Experiment - Comprehensive Report
==============================================================================
Generated: {timestamp}
VCF File: {self.vcf_file_path}

=== SYSTEM ENVIRONMENT ===
CPU Count: {psutil.cpu_count()}
Total Memory: {psutil.virtual_memory().total / (1024**3):.1f} GB
Python Version: {sys.version.split()[0]}
Platform: {sys.platform}

=== FILE ANALYSIS ===
"""
        
        if 'file_stats' in self.results:
            fs = self.results['file_stats']
            report += f"""File Size: {fs.get('file_size_mb', 0):.2f} MB
Total Lines: {fs.get('total_lines', 0):,}
Header Lines: {fs.get('header_lines', 0):,}
Variant Lines: {fs.get('variant_lines', 0):,}
"""
        
        report += f"""
=== DATA LOADING PERFORMANCE ===
Schema Creation Time: {self.results.get('schema_creation_time', 0):.2f} seconds
Data Loading Time: {self.results.get('load_time', 0):.2f} seconds
Sample Insert Time: {self.results.get('sample_insert_time', 0):.2f} seconds

=== DATA LOADED ===
Variants Loaded: {self.results.get('variant_count', 0):,}
Genotypes Loaded: {self.results.get('genotype_count', 0):,}
Samples: {self.results.get('sample_count', 0):,}
Loading Rate: {self.results.get('loading_rate_variants_per_sec', 0):.0f} variants/second
"""
        
        # Chromosome distribution
        if 'chromosome_distribution' in self.results:
            report += f"\nChromosome Distribution:\n"
            for chrom, count in sorted(self.results['chromosome_distribution'].items()):
                report += f"  {chrom}: {count:,} variants\n"
        
        # Quality statistics
        if 'quality_stats' in self.results:
            qs = self.results['quality_stats']
            report += f"""
Quality Statistics:
  Total with Quality Scores: {qs.get('count', 0):,}
  Mean Quality: {qs.get('mean', 0):.2f}
  Median Quality: {qs.get('median', 0):.2f}
  Min Quality: {qs.get('min', 0):.2f}
  Max Quality: {qs.get('max', 0):.2f}
"""
        
        # Genotype distribution
        if 'genotype_distribution' in self.results:
            report += f"\nGenotype Distribution:\n"
            top_genotypes = sorted(self.results['genotype_distribution'].items(), 
                                 key=lambda x: x[1], reverse=True)[:10]
            for gt, count in top_genotypes:
                report += f"  {gt}: {count:,}\n"
        
        # System resource usage during loading
        if 'system_stats_loading' in self.results:
            loading_stats = self.results['system_stats_loading']
            report += f"""
=== SYSTEM RESOURCES DURING LOADING ===
Duration: {loading_stats.get('duration_seconds', 0):.1f} seconds
CPU Usage:
  Average: {loading_stats.get('cpu_percent', {}).get('avg', 0):.1f}%
  Peak: {loading_stats.get('cpu_percent', {}).get('max', 0):.1f}%
Memory Usage:
  Average: {loading_stats.get('memory_percent', {}).get('avg', 0):.1f}%
  Peak: {loading_stats.get('memory_percent', {}).get('max', 0):.1f}%
  Process Memory Peak: {loading_stats.get('memory_used_gb', {}).get('max', 0):.2f} GB
Disk I/O:
  Data Read: {loading_stats.get('disk_read_mb', {}).get('final', 0):.1f} MB
  Data Written: {loading_stats.get('disk_write_mb', {}).get('final', 0):.1f} MB
Network I/O:
  Data Sent: {loading_stats.get('network_sent_mb', {}).get('final', 0):.1f} MB
  Data Received: {loading_stats.get('network_recv_mb', {}).get('final', 0):.1f} MB
"""
        
        # Database statistics
        if 'database_stats' in self.results:
            db_stats = self.results['database_stats']
            report += f"""
=== DATABASE STATISTICS ===
Total Database Size: {db_stats.get('database_size_mb', 0):.2f} MB
Storage Size: {db_stats.get('storage_size_mb', 0):.2f} MB
Index Size: {db_stats.get('index_size_mb', 0):.2f} MB
Total Objects: {db_stats.get('object_count', 0):,}

Collection Statistics:
"""
            for coll_name, coll_stats in db_stats.get('collection_stats', {}).items():
                if coll_stats:
                    report += f"""  {coll_name}:
    Document Count: {coll_stats.get('count', 0):,}
    Data Size: {coll_stats.get('size_mb', 0):.2f} MB
    Storage Size: {coll_stats.get('storage_size_mb', 0):.2f} MB
    Average Object Size: {coll_stats.get('avg_obj_size', 0):.0f} bytes
    Index Count: {coll_stats.get('index_count', 0)}
    Index Size: {coll_stats.get('total_index_size_mb', 0):.2f} MB
"""
        
        # Query performance results
        report += f"\n=== QUERY PERFORMANCE RESULTS ===\n"
        
        if 'queries' in self.results:
            for query_name, query_data in self.results['queries'].items():
                report += f"""
{query_data['description']}:
  Average Time: {query_data['avg_time']*1000:.2f} ms
  Min Time: {query_data['min_time']*1000:.2f} ms
  Max Time: {query_data['max_time']*1000:.2f} ms
  Std Deviation: {query_data['std_time']*1000:.2f} ms
  Results Returned: {query_data['result_count']:,}
  
  Sample Results: {str(query_data['sample_results'][:3])}
"""
        
        # System resources during queries
        if 'system_stats_queries' in self.results:
            query_stats = self.results['system_stats_queries']
            report += f"""
=== SYSTEM RESOURCES DURING QUERIES ===
Duration: {query_stats.get('duration_seconds', 0):.1f} seconds
CPU Usage:
  Average: {query_stats.get('cpu_percent', {}).get('avg', 0):.1f}%
  Peak: {query_stats.get('cpu_percent', {}).get('max', 0):.1f}%
Memory Usage:
  Average: {query_stats.get('memory_percent', {}).get('avg', 0):.1f}%
  Peak: {query_stats.get('memory_percent', {}).get('max', 0):.1f}%
"""
        
        # Performance summary
        if 'database_totals' in self.results:
            totals = self.results['database_totals']
            avg_query_time = 0
            if 'queries' in self.results and self.results['queries']:
                avg_query_time = statistics.mean([q['avg_time'] for q in self.results['queries'].values()])
            
            storage_efficiency = 0
            if 'database_stats' in self.results:
                db_size_mb = self.results['database_stats'].get('database_size_mb', 1)
                if db_size_mb > 0:
                    storage_efficiency = self.results.get('variant_count', 0) / db_size_mb
            
            report += f"""
=== PERFORMANCE SUMMARY ===
Total Data Points: {totals.get('total_variants', 0):,} variants, {totals.get('total_genotypes', 0):,} genotypes
Loading Efficiency: {self.results.get('loading_rate_variants_per_sec', 0):.0f} variants/second
Storage Efficiency: {storage_efficiency:.0f} variants/MB
Average Query Response: {avg_query_time * 1000:.1f} ms
"""
        
        report += f"""
=== MONGODB-SPECIFIC RECOMMENDATIONS ===
Based on the performance analysis:

1. Loading Performance:
   - Loading rate of {self.results.get('loading_rate_variants_per_sec', 0):.0f} variants/second
   - Consider using bulk operations for better throughput
   - MongoDB's document model handles VCF data naturally

2. Query Performance:"""
        
        if 'queries' in self.results and self.results['queries']:
            avg_query_time = statistics.mean([q['avg_time'] for q in self.results['queries'].values()])
            report += f"""
   - Average query time: {avg_query_time * 1000:.1f} ms
   - Consider compound indexes for range queries
   - Aggregation pipeline performs well for complex analyses"""
        
        report += f"""

3. Storage Optimization:
   - Database size: {self.results.get('database_stats', {}).get('database_size_mb', 0):.2f} MB
   - MongoDB's BSON format is efficient for nested VCF data
   - Consider sharding for very large datasets
   - Regular compaction recommended for write-heavy workloads

4. Schema Design:
   - Document model naturally fits VCF structure
   - Embedded sample data reduces joins
   - Consider separating large INFO fields if needed

==============================================================================
"""
        
        # Write to file
        with open(output_file, 'w') as f:
            f.write(report)
        
        # Also save as JSON for programmatic access
        json_file = output_file.replace('.txt', '.json')
        with open(json_file, 'w') as f:
            # Clean results for JSON serialization
            clean_results = clean_for_json(self.results)
            json.dump(clean_results, f, indent=2, default=str)
        
        print(f"✅ Comprehensive report saved to {output_file}")
        print(f"✅ JSON results saved to {json_file}")
        print(report)
        
        return report
    
    def create_enhanced_visualizations(self):
        """Create comprehensive performance visualization charts"""
        if 'queries' not in self.results or not self.results['queries']:
            print("No query results to visualize")
            return
        
        try:
            # Create a comprehensive dashboard
            fig = plt.figure(figsize=(20, 16))
            
            # 1. Query Performance (Top Left)
            ax1 = plt.subplot(3, 3, 1)
            queries = self.results['queries']
            query_names = list(queries.keys())
            avg_times = [queries[name]['avg_time'] * 1000 for name in query_names]
            
            bars = ax1.bar(range(len(query_names)), avg_times, alpha=0.7, color='steelblue')
            ax1.set_xlabel('Query Type')
            ax1.set_ylabel('Average Time (ms)')
            ax1.set_title('MongoDB Query Performance')
            ax1.set_xticks(range(len(query_names)))
            ax1.set_xticklabels([name.replace('_', '\n') for name in query_names], rotation=45, ha='right')
            if max(avg_times) > 0:
                ax1.set_yscale('log')
            
            # 2. System Resources During Loading (Top Center)
            ax2 = plt.subplot(3, 3, 2)
            if 'system_stats_loading' in self.results:
                loading_stats = self.results['system_stats_loading']
                
                metrics = ['CPU %', 'Memory %', 'Disk Read MB', 'Disk Write MB']
                values = [
                    loading_stats.get('cpu_percent', {}).get('avg', 0),
                    loading_stats.get('memory_percent', {}).get('avg', 0),
                    loading_stats.get('disk_read_mb', {}).get('final', 0),
                    loading_stats.get('disk_write_mb', {}).get('final', 0)
                ]
                
                bars = ax2.bar(metrics, values, alpha=0.7, color=['red', 'orange', 'green', 'blue'])
                ax2.set_ylabel('Value')
                ax2.set_title('System Resources During Loading')
                if max(values) > 0:
                    ax2.set_yscale('log')
                
                for bar, val in zip(bars, values):
                    height = bar.get_height()
                    ax2.text(bar.get_x() + bar.get_width()/2., height,
                            f'{val:.1f}', ha='center', va='bottom')
            else:
                ax2.text(0.5, 0.5, 'No loading stats available', ha='center', va='center', transform=ax2.transAxes)
                ax2.set_title('System Resources During Loading')
            
            # 3. Data Distribution (Top Right)
            ax3 = plt.subplot(3, 3, 3)
            if 'chromosome_distribution' in self.results and self.results['chromosome_distribution']:
                chrom_data = self.results['chromosome_distribution']
                chromosomes = list(chrom_data.keys())[:10]  # Top 10
                counts = [chrom_data[c] for c in chromosomes]
                
                ax3.bar(chromosomes, counts, alpha=0.7, color='green')
                ax3.set_xlabel('Chromosome')
                ax3.set_ylabel('Variant Count')
                ax3.set_title('Variant Distribution by Chromosome')
                plt.setp(ax3.get_xticklabels(), rotation=45)
            else:
                ax3.text(0.5, 0.5, 'No chromosome data available', ha='center', va='center', transform=ax3.transAxes)
                ax3.set_title('Variant Distribution by Chromosome')
            
            # 4. Database Size Breakdown (Middle Left)
            ax4 = plt.subplot(3, 3, 4)
            if ('database_stats' in self.results and 
                'collection_stats' in self.results['database_stats'] and 
                self.results['database_stats']['collection_stats']):
                
                coll_stats = self.results['database_stats']['collection_stats']
                coll_names = []
                coll_sizes = []
                
                for coll_name, stats in coll_stats.items():
                    if stats and 'size_mb' in stats:
                        coll_names.append(coll_name)
                        coll_sizes.append(stats['size_mb'])
                
                if coll_names and coll_sizes:
                    ax4.pie(coll_sizes, labels=coll_names, autopct='%1.1f%%', startangle=90)
                    ax4.set_title('Database Storage by Collection')
                else:
                    ax4.text(0.5, 0.5, 'No collection size data', ha='center', va='center', transform=ax4.transAxes)
                    ax4.set_title('Database Storage by Collection')
            else:
                ax4.text(0.5, 0.5, 'No database stats available', ha='center', va='center', transform=ax4.transAxes)
                ax4.set_title('Database Storage by Collection')
            
            # 5. Quality Distribution (Middle Center)
            ax5 = plt.subplot(3, 3, 5)
            if 'quality_stats' in self.results:
                qs = self.results['quality_stats']
                quality_metrics = ['Mean', 'Median', 'Min', 'Max']
                quality_values = [qs.get('mean', 0), qs.get('median', 0), qs.get('min', 0), qs.get('max', 0)]
                
                ax5.bar(quality_metrics, quality_values, alpha=0.7, color='purple')
                ax5.set_ylabel('Quality Score')
                ax5.set_title('Quality Score Statistics')
            else:
                ax5.text(0.5, 0.5, 'No quality data available', ha='center', va='center', transform=ax5.transAxes)
                ax5.set_title('Quality Score Statistics')
            
            # 6. Genotype Distribution (Middle Right)
            ax6 = plt.subplot(3, 3, 6)
            if 'genotype_distribution' in self.results and self.results['genotype_distribution']:
                gt_data = self.results['genotype_distribution']
                top_genotypes = sorted(gt_data.items(), key=lambda x: x[1], reverse=True)[:8]
                genotypes = [gt[0] for gt in top_genotypes]
                counts = [gt[1] for gt in top_genotypes]
                
                ax6.bar(genotypes, counts, alpha=0.7, color='orange')
                ax6.set_xlabel('Genotype')
                ax6.set_ylabel('Count')
                ax6.set_title('Top Genotype Distribution')
                plt.setp(ax6.get_xticklabels(), rotation=45)
            else:
                ax6.text(0.5, 0.5, 'No genotype data available', ha='center', va='center', transform=ax6.transAxes)
                ax6.set_title('Top Genotype Distribution')
            
            # 7. Loading Performance Timeline (Bottom Left)
            ax7 = plt.subplot(3, 3, 7)
            loading_metrics = ['Schema Creation', 'Data Loading', 'Sample Insert']
            loading_times = [
                self.results.get('schema_creation_time', 0),
                self.results.get('load_time', 0),
                self.results.get('sample_insert_time', 0)
            ]
            
            ax7.barh(loading_metrics, loading_times, alpha=0.7, color='cyan')
            ax7.set_xlabel('Time (seconds)')
            ax7.set_title('Loading Performance Breakdown')
            
            # 8. Query Performance Distribution (Bottom Center)
            ax8 = plt.subplot(3, 3, 8)
            if 'queries' in self.results and self.results['queries']:
                query_times = [q['avg_time'] * 1000 for q in self.results['queries'].values()]
                if query_times and max(query_times) > 0:
                    ax8.hist(query_times, bins=min(10, len(query_times)), alpha=0.7, color='lightblue', edgecolor='black')
                    ax8.set_xlabel('Query Time (ms)')
                    ax8.set_ylabel('Frequency')
                    ax8.set_title('Query Performance Distribution')
                    if max(query_times) > 10:
                        ax8.set_xscale('log')
                else:
                    ax8.text(0.5, 0.5, 'No query time data', ha='center', va='center', transform=ax8.transAxes)
                    ax8.set_title('Query Performance Distribution')
            else:
                ax8.text(0.5, 0.5, 'No query data available', ha='center', va='center', transform=ax8.transAxes)
                ax8.set_title('Query Performance Distribution')
            
            # 9. Summary Statistics (Bottom Right)
            ax9 = plt.subplot(3, 3, 9)
            ax9.axis('off')  # Turn off axis for text
            
            avg_query_time = 0
            if 'queries' in self.results and self.results['queries']:
                avg_query_time = statistics.mean([q['avg_time'] for q in self.results['queries'].values()])
            
            peak_cpu = 0
            peak_memory = 0
            if 'system_stats_loading' in self.results:
                peak_cpu = self.results['system_stats_loading'].get('cpu_percent', {}).get('max', 0)
                peak_memory = self.results['system_stats_loading'].get('memory_used_gb', {}).get('max', 0)
            
            db_size = 0
            if 'database_stats' in self.results:
                db_size = self.results['database_stats'].get('database_size_mb', 0)
            
            summary_text = f"""
MONGODB EXPERIMENT SUMMARY

Data Loaded:
• {self.results.get('variant_count', 0):,} variants
• {self.results.get('genotype_count', 0):,} genotypes  
• {self.results.get('sample_count', 0)} samples

Performance:
• Loading: {self.results.get('loading_rate_variants_per_sec', 0):.0f} var/sec
• Avg Query: {avg_query_time * 1000:.1f} ms
• DB Size: {db_size:.1f} MB

System Peak Usage:
• CPU: {peak_cpu:.1f}%
• Memory: {peak_memory:.1f} GB
"""
            
            ax9.text(0.1, 0.9, summary_text, transform=ax9.transAxes, fontsize=10,
                    verticalalignment='top', fontfamily='monospace',
                    bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
            
            plt.tight_layout()
            
            # Get file prefix for naming
            file_prefix = self._get_file_prefix()
            chart_filename = f'{file_prefix}_mongodb_vcf_comprehensive_performance.png'
            plt.savefig(chart_filename, dpi=300, bbox_inches='tight')
            plt.show()
            
            # Create a separate system resource timeline if we have the data
            if ('system_stats_loading' in self.results and 
                self.system_monitor.metrics['timestamps']):
                self.create_resource_timeline()
            
            print(f"✅ Enhanced visualizations saved as '{chart_filename}'")
            
        except Exception as e:
            print(f"Warning: Could not create visualizations: {e}")
            print("Continuing without visualizations...")
    
    def create_resource_timeline(self):
        """Create timeline charts for system resource usage"""
        try:
            metrics = self.system_monitor.metrics
            if not metrics['timestamps']:
                return
                
            # Convert timestamps to relative time
            start_time = metrics['timestamps'][0]
            relative_times = [(t - start_time) / 60 for t in metrics['timestamps']]  # Convert to minutes
            
            fig, axes = plt.subplots(2, 2, figsize=(15, 10))
            
            # CPU Usage
            if metrics['cpu_percent']:
                axes[0, 0].plot(relative_times, metrics['cpu_percent'], 'r-', linewidth=2)
            axes[0, 0].set_title('CPU Usage Over Time')
            axes[0, 0].set_xlabel('Time (minutes)')
            axes[0, 0].set_ylabel('CPU Usage (%)')
            axes[0, 0].grid(True, alpha=0.3)
            
            # Memory Usage
            if metrics['memory_used_gb']:
                axes[0, 1].plot(relative_times, metrics['memory_used_gb'], 'b-', linewidth=2)
            axes[0, 1].set_title('Process Memory Usage Over Time')
            axes[0, 1].set_xlabel('Time (minutes)')
            axes[0, 1].set_ylabel('Memory Usage (GB)')
            axes[0, 1].grid(True, alpha=0.3)
            
            # Disk I/O
            if metrics['disk_read_mb'] and metrics['disk_write_mb']:
                axes[1, 1].set_title('Network I/O Over Time')
            axes[1, 1].set_xlabel('Time (minutes)')
            axes[1, 1].set_ylabel('Network I/O (MB)')
            axes[1, 1].grid(True, alpha=0.3)
            
            plt.tight_layout()
            
            file_prefix = self._get_file_prefix()
            timeline_filename = f'{file_prefix}_mongodb_vcf_resource_timeline.png'
            plt.savefig(timeline_filename, dpi=300, bbox_inches='tight')
            plt.show()
            
            print(f"✅ Resource timeline saved as '{timeline_filename}'")
            
        except Exception as e:
            print(f"Warning: Could not create resource timeline: {e}")
    
    def _get_file_prefix(self):
        """Extract file prefix from VCF filename for output naming"""
        vcf_basename = os.path.basename(self.vcf_file_path)
        if vcf_basename.startswith('10k'):
            return '10k'
        elif vcf_basename.startswith('50k'):
            return '50k'
        elif vcf_basename.startswith('100k'):
            return '100k'
        else:
            # Default fallback
            return vcf_basename.replace('.vcf', '').replace('_test', '')
    
    def run_experiment(self):
        """Run the complete enhanced experiment"""
        file_prefix = self._get_file_prefix()
        print(f"=== Enhanced MongoDB VCF Performance Experiment ({file_prefix}) ===\n")
        
        # Test connection
        if not self.test_connection():
            return False
        
        try:
            # Create schema
            self.create_schema()
            
            # Load data with monitoring
            self.load_data()
            
            # Run performance tests
            self.run_performance_tests()
            
            # Get database statistics
            self.get_database_statistics()
            
            # Generate comprehensive report with appropriate filename
            report_filename = f"{file_prefix}_mongodb_vcf_comprehensive_report.txt"
            self.generate_comprehensive_report(report_filename)
            
            # Create enhanced visualizations
            self.create_enhanced_visualizations()
            
            print(f"\n✅ Enhanced experiment for {file_prefix} completed successfully!")
            return True
            
        except Exception as e:
            print(f"❌ Experiment failed: {e}")
            import traceback
            traceback.print_exc()
            return False
        finally:
            # Clean up MongoDB connection
            if self.client:
                self.client.close()

def run_multiple_experiments(vcf_files):
    """Run experiments on multiple VCF files"""
    all_results = {}
    
    for vcf_file in vcf_files:
        if not os.path.exists(vcf_file):
            print(f"Warning: VCF file '{vcf_file}' not found, skipping...")
            continue
            
        print(f"\n{'='*80}")
        print(f"STARTING EXPERIMENT FOR: {vcf_file}")
        print(f"{'='*80}")
        
        # Create experiment instance
        experiment = MongoDBVCFExperiment(vcf_file)
        
        # Run experiment
        success = experiment.run_experiment()
        
        if success:
            file_prefix = experiment._get_file_prefix()
            all_results[file_prefix] = experiment.results
            
            print(f"\n✅ Experiment for {vcf_file} completed successfully!")
            print(f"   Variants loaded: {experiment.results.get('variant_count', 0):,}")
            print(f"   Loading time: {experiment.results.get('load_time', 0):.2f} seconds")
            if experiment.results.get('queries'):
                avg_query_time = statistics.mean([q['avg_time'] for q in experiment.results['queries'].values()])
                print(f"   Average query time: {avg_query_time * 1000:.1f} ms")
            else:
                print("   No query results")
        else:
            print(f"❌ Experiment for {vcf_file} failed!")
    
    # Generate comparative report
    if len(all_results) > 1:
        generate_comparative_report(all_results)
    
    return all_results

def generate_comparative_report(all_results):
    """Generate a comparative report across all experiments"""
    print("\nGenerating comparative report...")
    
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    report = f"""
==============================================================================
MongoDB VCF Performance Comparison Report
==============================================================================
Generated: {timestamp}

=== COMPARATIVE ANALYSIS ===

"""
    
    # Create comparison table
    datasets = list(all_results.keys())
    
    report += f"{'Dataset':<10} {'Variants':<12} {'Genotypes':<12} {'Load Time':<12} {'Load Rate':<15} {'Avg Query':<12} {'DB Size':<12}\n"
    report += f"{'-'*10} {'-'*12} {'-'*12} {'-'*12} {'-'*15} {'-'*12} {'-'*12}\n"
    
    for dataset in sorted(datasets):
        results = all_results[dataset]
        variants = results.get('variant_count', 0)
        genotypes = results.get('genotype_count', 0)
        load_time = results.get('load_time', 0)
        load_rate = results.get('loading_rate_variants_per_sec', 0)
        
        avg_query_time = 0
        if 'queries' in results and results['queries']:
            avg_query_time = statistics.mean([q['avg_time'] for q in results['queries'].values()])
        
        db_size = results.get('database_stats', {}).get('database_size_mb', 0)
        
        report += f"{dataset:<10} {variants:<12,} {genotypes:<12,} {load_time:<12.1f} {load_rate:<15.0f} {avg_query_time*1000:<12.1f} {db_size:<12.1f}\n"
    
    # Performance scaling analysis
    report += f"""

=== SCALING ANALYSIS ===

Loading Performance Scaling:
"""
    
    for dataset in sorted(datasets):
        results = all_results[dataset]
        variants = results.get('variant_count', 0)
        load_time = results.get('load_time', 0)
        throughput = variants / load_time if load_time > 0 else 0
        
        report += f"  {dataset}: {throughput:.0f} variants/second\n"
    
    report += f"""

Query Performance Scaling:
"""
    
    for dataset in sorted(datasets):
        results = all_results[dataset]
        if 'queries' in results and results['queries']:
            avg_query_time = statistics.mean([q['avg_time'] for q in results['queries'].values()])
            report += f"  {dataset}: {avg_query_time*1000:.1f} ms average\n"
    
    report += f"""

Storage Efficiency:
"""
    
    for dataset in sorted(datasets):
        results = all_results[dataset]
        variants = results.get('variant_count', 0)
        db_size = results.get('database_stats', {}).get('database_size_mb', 1)
        efficiency = variants / db_size if db_size > 0 else 0
        
        report += f"  {dataset}: {efficiency:.0f} variants/MB\n"
    
    report += f"""

=== RECOMMENDATIONS ===

1. MongoDB Performance Characteristics:
   - Document model handles VCF data structure naturally
   - Loading performance scales reasonably with dataset size
   - Query performance varies by complexity and indexing

2. Optimization Strategies:
   - Use compound indexes for range queries
   - Consider sharding for datasets > 1M variants
   - Bulk operations provide better loading performance
   - Aggregation pipeline is efficient for complex analyses

3. Scaling Observations:
   - Memory usage grows with dataset size
   - Index size becomes significant factor
   - Network I/O minimal due to local storage

==============================================================================
"""
    
    # Write comparative report
    with open('mongodb_vcf_comparative_report.txt', 'w') as f:
        f.write(report)
    
    # Save comparative data as JSON
    with open('mongodb_vcf_comparative_results.json', 'w') as f:
        clean_results = clean_for_json(all_results)
        json.dump(clean_results, f, indent=2, default=str)
    
    print("✅ Comparative report saved to 'mongodb_vcf_comparative_report.txt'")
    print("✅ Comparative data saved to 'mongodb_vcf_comparative_results.json'")
    print(report)

def main():
    """Main function"""
    # Define the VCF files to test
    vcf_files = ['10k_test.vcf', '50k_test.vcf', '100k_test.vcf']
    
    print("MongoDB VCF Performance Experiment Suite")
    print("=========================================")
    print(f"Testing files: {', '.join(vcf_files)}")
    
    # Check if files exist
    missing_files = [f for f in vcf_files if not os.path.exists(f)]
    if missing_files:
        print(f"\nError: The following VCF files were not found:")
        for f in missing_files:
            print(f"  - {f}")
        print("\nPlease ensure all test files are in the current directory.")
        sys.exit(1)
    
    # Run experiments on all files
    all_results = run_multiple_experiments(vcf_files)
    
    if all_results:
        print(f"\n{'='*80}")
        print("ALL EXPERIMENTS COMPLETED")
        print(f"{'='*80}")
        print(f"Successfully processed {len(all_results)} datasets")
        
        for dataset, results in all_results.items():
            print(f"\n{dataset} Results:")
            print(f"  - Variants: {results.get('variant_count', 0):,}")
            print(f"  - Load time: {results.get('load_time', 0):.2f} seconds")
            print(f"  - DB size: {results.get('database_stats', {}).get('database_size_mb', 0):.1f} MB")
        
        print(f"\nGenerated files:")
        for dataset in all_results.keys():
            print(f"  - {dataset}_mongodb_vcf_comprehensive_report.txt")
            print(f"  - {dataset}_mongodb_vcf_comprehensive_report.json") 
            print(f"  - {dataset}_mongodb_vcf_comprehensive_performance.png")
            print(f"  - {dataset}_mongodb_vcf_resource_timeline.png")
        
        print(f"  - mongodb_vcf_comparative_report.txt")
        print(f"  - mongodb_vcf_comparative_results.json")
        
    else:
        print("❌ No experiments completed successfully")
        sys.exit(1)

if __name__ == "__main__":
    main()
