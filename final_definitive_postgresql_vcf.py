#!/usr/bin/env python3
"""
Enhanced PostgreSQL VCF Performance Experiment
Standalone script to load VCF data into PostgreSQL and measure query performance
with comprehensive system resource monitoring
"""

import psycopg2
from psycopg2.extras import Json
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

class PostgreSQLVCFExperiment:
    """Main experiment class for PostgreSQL VCF performance testing"""
    
    def __init__(self, vcf_file_path, db_params=None):
        self.vcf_file_path = vcf_file_path
        self.db_params = db_params or {
            'host': 'localhost',
            'database': 'vcf_db',
            'user': 'postgres',
            'password': 'password',
            'port': 5432
        }
        self.results = {}
        self.system_monitor = SystemMonitor()
        
    def test_connection(self):
        """Test database connection"""
        try:
            conn = psycopg2.connect(**self.db_params)
            cur = conn.cursor()
            cur.execute("SELECT version();")
            version = cur.fetchone()[0]
            print(f"✅ PostgreSQL connected successfully!")
            print(f"   Version: {version}")
            cur.close()
            conn.close()
            return True
        except Exception as e:
            print(f"❌ PostgreSQL connection failed: {e}")
            print("Make sure PostgreSQL is running:")
            print("docker run --name postgres-vcf -e POSTGRES_PASSWORD=password -e POSTGRES_DB=vcf_db -p 5432:5432 -d postgres:15")
            return False
    
    def create_schema(self):
        """Create database schema"""
        print("Creating PostgreSQL schema...")
        
        start_time = time.time()
        
        conn = psycopg2.connect(**self.db_params)
        cur = conn.cursor()
        
        try:
            # Drop existing tables
            cur.execute("DROP TABLE IF EXISTS genotypes CASCADE;")
            cur.execute("DROP TABLE IF EXISTS variants CASCADE;")
            cur.execute("DROP TABLE IF EXISTS samples CASCADE;")
            
            # Create tables
            schema_sql = """
            CREATE TABLE variants (
                id SERIAL PRIMARY KEY,
                chrom VARCHAR(1000) NOT NULL,
                pos BIGINT NOT NULL,
                variant_id VARCHAR(255),
                ref TEXT NOT NULL,
                alt TEXT NOT NULL,
                qual FLOAT,
                filter VARCHAR(255),
                info JSONB,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            );

            CREATE TABLE samples (
                id SERIAL PRIMARY KEY,
                sample_name VARCHAR(255) UNIQUE NOT NULL
            );

            CREATE TABLE genotypes (
                id SERIAL PRIMARY KEY,
                variant_id INTEGER REFERENCES variants(id),
                sample_id INTEGER REFERENCES samples(id),
                genotype VARCHAR(10),
                genotype_quality INTEGER,
                depth INTEGER,
                other_format JSONB
            );
            """
            
            cur.execute(schema_sql)
            conn.commit()
            
            # Verify tables were created
            cur.execute("SELECT EXISTS (SELECT FROM pg_tables WHERE schemaname = 'public' AND tablename = 'variants');")
            variants_exists = cur.fetchone()[0]
            if not variants_exists:
                raise RuntimeError("Failed to create variants table")
            print("✅ Variants table verified in database")
            
            # Create indexes
            index_sql = """
            CREATE INDEX idx_variants_chrom_pos ON variants(chrom, pos);
            CREATE INDEX idx_variants_qual ON variants(qual);
            CREATE INDEX idx_variants_chrom ON variants(chrom);
            CREATE INDEX idx_genotypes_variant_sample ON genotypes(variant_id, sample_id);
            CREATE INDEX idx_info_gin ON variants USING GIN (info);
            """
            
            cur.execute(index_sql)
            conn.commit()
            
        finally:
            cur.close()
            conn.close()
        
        schema_time = time.time() - start_time
        self.results['schema_creation_time'] = schema_time
        
        print(f"✅ Schema created successfully! ({schema_time:.2f}s)")
    
    def load_data(self):
        """Load VCF data into PostgreSQL with resource monitoring"""
        print(f"Loading VCF data from {self.vcf_file_path}...")
        
        if not os.path.exists(self.vcf_file_path):
            raise FileNotFoundError(f"VCF file not found: {self.vcf_file_path}")
        
        # Analyze file first
        vcf_parser = VCFParser(self.vcf_file_path)
        file_stats = vcf_parser.analyze_file()
        self.results['file_stats'] = file_stats
        
        # Start system monitoring
        self.system_monitor.start_monitoring()
        
        conn = psycopg2.connect(**self.db_params)
        cur = conn.cursor()
        
        try:
            # Parse VCF file
            vcf_parser.parse_header()
            
            print(f"Found {len(vcf_parser.samples)} samples in VCF file")
            
            # Insert samples
            sample_insert_start = time.time()
            for sample in vcf_parser.samples:
                cur.execute(
                    "INSERT INTO samples (sample_name) VALUES (%s) ON CONFLICT (sample_name) DO NOTHING",
                    (sample,)
                )
            
            conn.commit()
            sample_insert_time = time.time() - sample_insert_start
            
            # Get sample IDs
            cur.execute("SELECT id, sample_name FROM samples")
            sample_id_map = {name: id for id, name in cur.fetchall()}
            
            # Track loading progress
            start_time = time.time()
            variant_count = 0
            genotype_count = 0
            batch_size = 1000
            progress_interval = 5000
            
            # Statistics tracking
            chrom_stats = defaultdict(int)
            qual_distribution = []
            genotype_stats = defaultdict(int)
            
            variant_batch = []
            
            # Store variants first
            for variant in vcf_parser.parse_variants():
                if variant['pos'] is None:
                    continue
                
                # Collect statistics
                chrom_stats[variant['chrom']] += 1
                if variant['qual'] is not None:
                    qual_distribution.append(variant['qual'])
                
                # Collect genotype statistics
                for sample_name, sample_data in variant['samples'].items():
                    if sample_name in sample_id_map:
                        gt = sample_data.get('GT', 'Unknown')
                        genotype_stats[gt] += 1
                
                # Prepare variant for batch insert
                variant_data = (
                    variant['chrom'],
                    variant['pos'],
                    variant['id'],
                    variant['ref'],
                    ','.join(variant['alt']),
                    variant['qual'],
                    ','.join(variant['filter']) if variant['filter'] else None,
                    Json(variant['info'])
                )
                variant_batch.append(variant_data)
                
                variant_count += 1
                
                # Batch insert when batch is full
                if len(variant_batch) >= batch_size:
                    # Insert variants
                    cur.executemany(
                        """INSERT INTO variants (chrom, pos, variant_id, ref, alt, qual, filter, info)
                           VALUES (%s, %s, %s, %s, %s, %s, %s, %s)""",
                        variant_batch
                    )
                    conn.commit()
                    variant_batch = []
                    
                    if variant_count % progress_interval == 0:
                        current_time = time.time()
                        elapsed = current_time - start_time
                        rate = variant_count / elapsed
                        print(f"   Progress: {variant_count:,} variants loaded ({rate:.0f} variants/sec)")
            
            # Insert remaining variants
            if variant_batch:
                cur.executemany(
                    """INSERT INTO variants (chrom, pos, variant_id, ref, alt, qual, filter, info)
                       VALUES (%s, %s, %s, %s, %s, %s, %s, %s)""",
                    variant_batch
                )
                conn.commit()
            
            # Now process genotypes in a separate pass
            print("   Inserting genotypes...")
            genotype_start_time = time.time()
            
            genotype_batch = []
            batch_count = 0
            
            # Re-parse file for genotypes
            for variant in vcf_parser.parse_variants():
                if variant['pos'] is None:
                    continue
                    
                # Get variant ID - use a more specific query to handle duplicates
                cur.execute(
                    """SELECT id FROM variants 
                       WHERE chrom = %s AND pos = %s AND ref = %s AND alt = %s 
                       LIMIT 1""",
                    (variant['chrom'], variant['pos'], variant['ref'], ','.join(variant['alt']))
                )
                result = cur.fetchone()
                if not result:
                    continue
                variant_id = result[0]
                
                # Prepare genotype batch
                for sample_name, sample_data in variant['samples'].items():
                    if sample_name in sample_id_map:
                        gq = sample_data.get('GQ', '')
                        dp = sample_data.get('DP', '')
                        
                        # Handle numeric conversion safely
                        gq_val = None
                        if gq and gq != '.' and gq.replace('.', '').isdigit():
                            try:
                                gq_val = int(float(gq))
                            except ValueError:
                                gq_val = None
                                
                        dp_val = None
                        if dp and dp != '.' and dp.replace('.', '').isdigit():
                            try:
                                dp_val = int(float(dp))
                            except ValueError:
                                dp_val = None
                        
                        genotype_data = (
                            variant_id,
                            sample_id_map[sample_name],
                            sample_data.get('GT'),
                            gq_val,
                            dp_val,
                            Json(sample_data)
                        )
                        genotype_batch.append(genotype_data)
                        genotype_count += 1
                
                # Insert genotype batch
                if len(genotype_batch) >= batch_size:
                    cur.executemany(
                        """INSERT INTO genotypes (variant_id, sample_id, genotype, genotype_quality, depth, other_format)
                           VALUES (%s, %s, %s, %s, %s, %s)""",
                        genotype_batch
                    )
                    conn.commit()
                    genotype_batch = []
                    batch_count += 1
                    
                    if batch_count % 10 == 0:
                        print(f"   Genotypes progress: {genotype_count:,}")
            
            # Insert remaining genotypes
            if genotype_batch:
                cur.executemany(
                    """INSERT INTO genotypes (variant_id, sample_id, genotype, genotype_quality, depth, other_format)
                       VALUES (%s, %s, %s, %s, %s, %s)""",
                    genotype_batch
                )
                conn.commit()
            
        finally:
            cur.close()
            conn.close()
        
        # Stop monitoring and collect stats
        self.system_monitor.stop_monitoring()
        load_time = time.time() - start_time
        genotype_time = time.time() - genotype_start_time
        
        # Store comprehensive results
        self.results.update({
            'load_time': load_time,
            'sample_insert_time': sample_insert_time,
            'genotype_insert_time': genotype_time,
            'variant_count': variant_count,
            'genotype_count': genotype_count,
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
        print(f"   Loaded {genotype_count:,} genotypes")
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
        
        conn = psycopg2.connect(**self.db_params)
        cur = conn.cursor()
        
        try:
            # Get database statistics first
            cur.execute("SELECT COUNT(*) FROM variants")
            total_variants = cur.fetchone()[0]
            
            cur.execute("SELECT COUNT(*) FROM genotypes")
            total_genotypes = cur.fetchone()[0]
            
            cur.execute("SELECT COUNT(*) FROM samples")
            total_samples = cur.fetchone()[0]
            
            # Get some sample data for parameterized queries
            cur.execute("SELECT DISTINCT chrom FROM variants LIMIT 1")
            sample_chrom = cur.fetchone()[0] if cur.rowcount > 0 else '1'
            
            cur.execute("SELECT sample_name FROM samples LIMIT 1")
            sample_name = cur.fetchone()[0] if cur.rowcount > 0 else None
            
            cur.execute("SELECT MIN(pos), MAX(pos) FROM variants WHERE chrom = %s", (sample_chrom,))
            pos_range = cur.fetchone()
            min_pos = pos_range[0] if pos_range and pos_range[0] else 1000000
            max_pos = pos_range[1] if pos_range and pos_range[1] else 2000000
            
            queries = {
                "total_variants": {
                    "sql": "SELECT COUNT(*) FROM variants",
                    "params": None,
                    "description": "Count total variants",
                    "expected_result_type": "count"
                },
                "total_genotypes": {
                    "sql": "SELECT COUNT(*) FROM genotypes", 
                    "params": None,
                    "description": "Count total genotypes",
                    "expected_result_type": "count"
                },
                "range_query": {
                    "sql": "SELECT COUNT(*) FROM variants WHERE chrom = %s AND pos BETWEEN %s AND %s",
                    "params": (sample_chrom, min_pos, min_pos + 100000),
                    "description": f"Count variants in range {sample_chrom}:{min_pos}-{min_pos + 100000}",
                    "expected_result_type": "count"
                },
                "quality_filter": {
                    "sql": "SELECT COUNT(*) FROM variants WHERE qual > %s",
                    "params": (30,),
                    "description": "Count high-quality variants (QUAL > 30)",
                    "expected_result_type": "count"
                },
                "quality_stats": {
                    "sql": "SELECT AVG(qual), MIN(qual), MAX(qual), STDDEV(qual) FROM variants WHERE qual IS NOT NULL",
                    "params": None,
                    "description": "Quality statistics",
                    "expected_result_type": "stats"
                },
                "chromosome_summary": {
                    "sql": "SELECT chrom, COUNT(*) as variant_count, AVG(qual) as avg_qual FROM variants GROUP BY chrom ORDER BY chrom",
                    "params": None,
                    "description": "Variants per chromosome with average quality",
                    "expected_result_type": "summary"
                },
                "info_field_query": {
                    "sql": "SELECT COUNT(*) FROM variants WHERE info ? 'DP'",
                    "params": None,
                    "description": "Count variants with depth information",
                    "expected_result_type": "count"
                },
                "genotype_distribution": {
                    "sql": "SELECT genotype, COUNT(*) FROM genotypes WHERE genotype IS NOT NULL GROUP BY genotype ORDER BY COUNT(*) DESC",
                    "params": None,
                    "description": "Genotype distribution",
                    "expected_result_type": "distribution"
                },
                "complex_join": {
                    "sql": """
                        SELECT v.chrom, v.pos, COUNT(g.id) as sample_count, 
                               AVG(g.depth) as avg_depth
                        FROM variants v 
                        JOIN genotypes g ON v.id = g.variant_id 
                        WHERE v.chrom = %s 
                        GROUP BY v.chrom, v.pos 
                        ORDER BY v.pos 
                        LIMIT 100
                    """,
                    "params": (sample_chrom,),
                    "description": f"Complex join query for chromosome {sample_chrom}",
                    "expected_result_type": "complex"
                },
                "depth_statistics": {
                    "sql": "SELECT AVG(depth), MIN(depth), MAX(depth), COUNT(*) FROM genotypes WHERE depth IS NOT NULL",
                    "params": None,
                    "description": "Depth statistics across all genotypes",
                    "expected_result_type": "stats"
                }
            }
            
            # Add sample-specific query if we have samples
            if sample_name:
                queries["sample_genotypes"] = {
                    "sql": """
                        SELECT v.chrom, v.pos, g.genotype, g.genotype_quality
                        FROM variants v 
                        JOIN genotypes g ON v.id = g.variant_id 
                        JOIN samples s ON g.sample_id = s.id 
                        WHERE s.sample_name = %s 
                        LIMIT 100
                    """,
                    "params": (sample_name,),
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
                        if query_info['params']:
                            cur.execute(query_info['sql'], query_info['params'])
                        else:
                            cur.execute(query_info['sql'])
                        
                        result = cur.fetchall()
                        query_time = time.time() - start_time
                        times.append(query_time)
                        
                        if run == 0:  # Store result from first run
                            result_count = len(result)
                            result_data = result[:10] if result else []
                            
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
                        'sample_results': result_data,
                        'times': valid_times,
                        'query_type': query_info['expected_result_type']
                    }
                    print(f"    Average time: {statistics.mean(valid_times)*1000:.2f} ms")
                    print(f"    Results returned: {result_count:,}")
                else:
                    print(f"    Query failed on all attempts")
                    
        finally:
            cur.close()
            conn.close()
        
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
        
        conn = psycopg2.connect(**self.db_params)
        cur = conn.cursor()
        
        try:
            stats = {}
            
            # Database size
            cur.execute("SELECT pg_size_pretty(pg_database_size(%s))", (self.db_params['database'],))
            db_size = cur.fetchone()[0]
            stats['database_size'] = db_size
            
            # Get actual size in bytes for calculations
            cur.execute("SELECT pg_database_size(%s)", (self.db_params['database'],))
            db_size_bytes = cur.fetchone()[0]
            stats['database_size_bytes'] = db_size_bytes
            stats['database_size_mb'] = db_size_bytes / (1024 * 1024)
            
            # Table sizes with row counts - use separate connection for each query block
            try:
                cur.execute("""
                    SELECT 
                        t.relname as table_name,
                        pg_size_pretty(pg_total_relation_size(t.oid)) as total_size,
                        pg_size_pretty(pg_relation_size(t.oid)) as table_size,
                        pg_size_pretty(pg_total_relation_size(t.oid) - pg_relation_size(t.oid)) as index_size,
                        pg_total_relation_size(t.oid) as total_size_bytes,
                        COALESCE(s.n_tup_ins, 0) as inserts,
                        COALESCE(s.n_tup_upd, 0) as updates,
                        COALESCE(s.n_tup_del, 0) as deletes,
                        COALESCE(s.n_live_tup, 0) as live_tuples,
                        COALESCE(s.n_dead_tup, 0) as dead_tuples
                    FROM pg_class t
                    LEFT JOIN pg_stat_user_tables s ON t.relname = s.relname
                    WHERE t.relkind = 'r' AND t.relname IN ('variants', 'genotypes', 'samples')
                    ORDER BY pg_total_relation_size(t.oid) DESC
                """)
                
                table_stats = cur.fetchall()
                stats['table_stats'] = table_stats
            except Exception as e:
                print(f"Warning: Could not get table stats: {e}")
                stats['table_stats'] = []
                conn.rollback()  # Roll back failed transaction
            
            # Index usage and performance - fresh transaction
            try:
                cur.execute("""
                    SELECT 
                        indexrelname as index_name,
                        relname as table_name,
                        COALESCE(idx_tup_read, 0) as idx_tup_read,
                        COALESCE(idx_tup_fetch, 0) as idx_tup_fetch,
                        pg_size_pretty(pg_relation_size(indexrelid)) as index_size,
                        pg_relation_size(indexrelid) as index_size_bytes
                    FROM pg_stat_user_indexes 
                    WHERE relname IN ('variants', 'genotypes', 'samples')
                    ORDER BY idx_tup_read DESC
                """)
                
                index_stats = cur.fetchall()
                stats['index_stats'] = index_stats
            except Exception as e:
                print(f"Warning: Could not get index stats: {e}")
                stats['index_stats'] = []
                conn.rollback()  # Roll back failed transaction
            
            # Query performance statistics (optional - requires pg_stat_statements)
            try:
                cur.execute("""
                    SELECT 
                        query,
                        calls,
                        total_time,
                        mean_time,
                        rows
                    FROM pg_stat_statements 
                    WHERE query LIKE '%variants%' OR query LIKE '%genotypes%' OR query LIKE '%samples%'
                    ORDER BY total_time DESC
                    LIMIT 10
                """) 
                query_stats = cur.fetchall()
                stats['query_performance'] = query_stats
            except Exception as e:
                # pg_stat_statements might not be enabled
                stats['query_performance'] = []
                conn.rollback()  # Roll back failed transaction
            
            # Connection and activity stats - fresh transaction  
            try:
                cur.execute("""
                    SELECT 
                        datname,
                        COALESCE(numbackends, 0),
                        COALESCE(xact_commit, 0),
                        COALESCE(xact_rollback, 0),
                        COALESCE(blks_read, 0),
                        COALESCE(blks_hit, 0),
                        COALESCE(tup_returned, 0),
                        COALESCE(tup_fetched, 0),
                        COALESCE(tup_inserted, 0),
                        COALESCE(tup_updated, 0),
                        COALESCE(tup_deleted, 0)
                    FROM pg_stat_database 
                    WHERE datname = %s
                """, (self.db_params['database'],))
                
                db_activity = cur.fetchone()
                if db_activity:
                    stats['database_activity'] = {
                        'connections': db_activity[1],
                        'commits': db_activity[2],
                        'rollbacks': db_activity[3],
                        'blocks_read': db_activity[4],
                        'blocks_hit': db_activity[5],
                        'cache_hit_ratio': (db_activity[5] / max(db_activity[4] + db_activity[5], 1)) * 100,
                        'tuples_returned': db_activity[6],
                        'tuples_fetched': db_activity[7],
                        'tuples_inserted': db_activity[8],
                        'tuples_updated': db_activity[9],
                        'tuples_deleted': db_activity[10]
                    }
            except Exception as e:
                print(f"Warning: Could not get database activity stats: {e}")
                stats['database_activity'] = {}
                conn.rollback()  # Roll back failed transaction
                
        finally:
            cur.close()
            conn.close()
        
        self.results['database_stats'] = stats
        print(f"✅ Database size: {stats.get('database_size', 'Unknown')}")
        print(f"   Cache hit ratio: {stats.get('database_activity', {}).get('cache_hit_ratio', 0):.1f}%")
        
        return stats
    
    def generate_comprehensive_report(self, output_file="postgresql_vcf_comprehensive_report.txt"):
        """Generate a detailed comprehensive report of the experiment"""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        report = f"""
==============================================================================
PostgreSQL VCF Performance Experiment - Comprehensive Report
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
Genotype Insert Time: {self.results.get('genotype_insert_time', 0):.2f} seconds

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
Total Database Size: {db_stats.get('database_size', 'Unknown')}
Database Size (MB): {db_stats.get('database_size_mb', 0):.2f}

Table Statistics:
"""
            for table_data in db_stats.get('table_stats', []):
                table_name, total_size, table_size, index_size, total_bytes, inserts, updates, deletes, live_tuples, dead_tuples = table_data
                report += f"""  {table_name}:
    Total Size: {total_size}
    Table Size: {table_size}
    Index Size: {index_size}
    Live Tuples: {live_tuples:,}
    Dead Tuples: {dead_tuples:,}
    Inserts: {inserts:,}
"""
            
            if 'database_activity' in db_stats:
                activity = db_stats['database_activity']
                report += f"""
Database Activity:
  Cache Hit Ratio: {activity.get('cache_hit_ratio', 0):.2f}%
  Total Commits: {activity.get('commits', 0):,}
  Total Rollbacks: {activity.get('rollbacks', 0):,}
  Blocks Read: {activity.get('blocks_read', 0):,}
  Blocks Hit: {activity.get('blocks_hit', 0):,}
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
=== RECOMMENDATIONS ===
Based on the performance analysis:

1. Loading Performance:
   - Loading rate of {self.results.get('loading_rate_variants_per_sec', 0):.0f} variants/second
   - Consider batch size optimization for better throughput
   - Monitor memory usage during large file processing

2. Query Performance:"""
        
        if 'queries' in self.results and self.results['queries']:
            avg_query_time = statistics.mean([q['avg_time'] for q in self.results['queries'].values()])
            report += f"""
   - Average query time: {avg_query_time * 1000:.1f} ms
   - Consider additional indexes for frequently queried columns"""
        
        cache_hit_ratio = 0
        if 'database_stats' in self.results:
            cache_hit_ratio = self.results['database_stats'].get('database_activity', {}).get('cache_hit_ratio', 0)
        
        report += f"""
   - Monitor cache hit ratio: {cache_hit_ratio:.1f}%

3. Storage Optimization:
   - Database size: {self.results.get('database_stats', {}).get('database_size', 'Unknown')}
   - Consider compression for large INFO fields
   - Regular VACUUM operations recommended

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
            ax1.set_title('Query Performance')
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
                'table_stats' in self.results['database_stats'] and 
                self.results['database_stats']['table_stats']):
                table_stats = self.results['database_stats']['table_stats']
                table_names = [t[0] for t in table_stats]
                table_sizes = [t[4] / (1024*1024) for t in table_stats]  # Convert to MB
                
                ax4.pie(table_sizes, labels=table_names, autopct='%1.1f%%', startangle=90)
                ax4.set_title('Database Storage by Table')
            else:
                ax4.text(0.5, 0.5, 'No database stats available', ha='center', va='center', transform=ax4.transAxes)
                ax4.set_title('Database Storage by Table')
            
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
            loading_metrics = ['Schema Creation', 'Data Loading', 'Sample Insert', 'Genotype Insert']
            loading_times = [
                self.results.get('schema_creation_time', 0),
                self.results.get('load_time', 0),
                self.results.get('sample_insert_time', 0),
                self.results.get('genotype_insert_time', 0)
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
            
            summary_text = f"""
EXPERIMENT SUMMARY

Data Loaded:
• {self.results.get('variant_count', 0):,} variants
• {self.results.get('genotype_count', 0):,} genotypes  
• {self.results.get('sample_count', 0)} samples

Performance:
• Loading: {self.results.get('loading_rate_variants_per_sec', 0):.0f} var/sec
• Avg Query: {avg_query_time * 1000:.1f} ms
• DB Size: {self.results.get('database_stats', {}).get('database_size', 'Unknown')}

System Peak Usage:
• CPU: {peak_cpu:.1f}%
• Memory: {peak_memory:.1f} GB
"""
            
            ax9.text(0.1, 0.9, summary_text, transform=ax9.transAxes, fontsize=10,
                    verticalalignment='top', fontfamily='monospace',
                    bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
            
            plt.tight_layout()
            plt.savefig('postgresql_vcf_comprehensive_performance.png', dpi=300, bbox_inches='tight')
            plt.show()
            
            # Create a separate system resource timeline if we have the data
            if ('system_stats_loading' in self.results and 
                self.system_monitor.metrics['timestamps']):
                self.create_resource_timeline()
            
            print("✅ Enhanced visualizations saved as 'postgresql_vcf_comprehensive_performance.png'")
            
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
                axes[1, 0].plot(relative_times, metrics['disk_read_mb'], 'g-', linewidth=2, label='Read')
                axes[1, 0].plot(relative_times, metrics['disk_write_mb'], 'orange', linewidth=2, label='Write')
                axes[1, 0].legend()
            axes[1, 0].set_title('Disk I/O Over Time')
            axes[1, 0].set_xlabel('Time (minutes)')
            axes[1, 0].set_ylabel('Disk I/O (MB)')
            axes[1, 0].grid(True, alpha=0.3)
            
            # Network I/O
            if metrics['network_sent_mb'] and metrics['network_recv_mb']:
                axes[1, 1].plot(relative_times, metrics['network_sent_mb'], 'purple', linewidth=2, label='Sent')
                axes[1, 1].plot(relative_times, metrics['network_recv_mb'], 'brown', linewidth=2, label='Received')
                axes[1, 1].legend()
            axes[1, 1].set_title('Network I/O Over Time')
            axes[1, 1].set_xlabel('Time (minutes)')
            axes[1, 1].set_ylabel('Network I/O (MB)')
            axes[1, 1].grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig('postgresql_vcf_resource_timeline.png', dpi=300, bbox_inches='tight')
            plt.show()
            
            print("✅ Resource timeline saved as 'postgresql_vcf_resource_timeline.png'")
            
        except Exception as e:
            print(f"Warning: Could not create resource timeline: {e}")
    
    def run_experiment(self):
        """Run the complete enhanced experiment"""
        print("=== Enhanced PostgreSQL VCF Performance Experiment ===\n")
        
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
            
            # Generate comprehensive report
            self.generate_comprehensive_report()
            
            # Create enhanced visualizations
            self.create_enhanced_visualizations()
            
            print("\n✅ Enhanced experiment completed successfully!")
            return True
            
        except Exception as e:
            print(f"❌ Experiment failed: {e}")
            import traceback
            traceback.print_exc()
            return False

def main():
    """Main function"""
    if len(sys.argv) != 2:
        print("Usage: python enhanced_postgresql_vcf_experiment.py <vcf_file_path>")
        print("Example: python enhanced_postgresql_vcf_experiment.py test.vcf")
        sys.exit(1)
    
    vcf_file = sys.argv[1]
    
    if not os.path.exists(vcf_file):
        print(f"Error: VCF file '{vcf_file}' not found")
        sys.exit(1)
    
    # Create and run experiment
    experiment = PostgreSQLVCFExperiment(vcf_file)
    success = experiment.run_experiment()
    
    if success:
        print("\n" + "="*80)
        print("EXPERIMENT RESULTS SUMMARY")
        print("="*80)
        print(f"VCF file processed: {vcf_file}")
        print(f"File size: {experiment.results.get('file_stats', {}).get('file_size_mb', 0):.2f} MB")
        print(f"Variants loaded: {experiment.results.get('variant_count', 0):,}")
        print(f"Genotypes loaded: {experiment.results.get('genotype_count', 0):,}")
        print(f"Samples processed: {experiment.results.get('sample_count', 0)}")
        print(f"Loading time: {experiment.results.get('load_time', 0):.2f} seconds")
        print(f"Loading rate: {experiment.results.get('loading_rate_variants_per_sec', 0):.0f} variants/second")
        print(f"Database size: {experiment.results.get('database_stats', {}).get('database_size', 'Unknown')}")
        
        # Calculate average query time safely
        avg_query_time = 0
        if 'queries' in experiment.results and experiment.results['queries']:
            avg_query_time = statistics.mean([q['avg_time'] for q in experiment.results['queries'].values()])
        print(f"Average query time: {avg_query_time * 1000:.1f} ms")
        
        # System resource summary
        if 'system_stats_loading' in experiment.results:
            loading_stats = experiment.results['system_stats_loading']
            print(f"Peak CPU usage: {loading_stats.get('cpu_percent', {}).get('max', 0):.1f}%")
            print(f"Peak memory usage: {loading_stats.get('memory_used_gb', {}).get('max', 0):.2f} GB")
        
        print(f"\nReports generated:")
        print(f"- Comprehensive report: postgresql_vcf_comprehensive_report.txt")
        print(f"- JSON results: postgresql_vcf_comprehensive_report.json")
        print(f"- Performance charts: postgresql_vcf_comprehensive_performance.png")
        print(f"- Resource timeline: postgresql_vcf_resource_timeline.png")
        print("="*80)
    else:
        print("Experiment failed. Check the error messages above.")
        sys.exit(1)

if __name__ == "__main__":
    main()