#!/usr/bin/env python3
"""
Test script for auto-queue selection functionality
"""

import tempfile
import os
import sys

# Add the parent directory to Python path
sys.path.insert(0, '/g/data/a56/software/qsub_tools')

def test_auto_queue_selection():
    """Test the auto-queue selection with different resource requirements"""
    
    # Create a minimal platform definition for testing
    platform_config = """
name: test_platform
type: pbs
host: test.example.com
description: Test platform for auto-queue selection
queues:
  - name: small
    type: compute
    description: Small jobs
    limits:
      walltime: 1h
      memory: 4GB
      ncpus: 4
    cost_per_hour: 1.0
    
  - name: normal
    type: compute
    description: Normal jobs  
    limits:
      walltime: 48h
      memory: 192GB
      ncpus: 48
    cost_per_hour: 2.0
    
  - name: hugemem
    type: compute
    description: High memory jobs
    limits:
      walltime: 48h
      memory: 1024GB
      ncpus: 64
    cost_per_hour: 8.0
    
  - name: gpu
    type: gpu
    description: GPU jobs
    limits:
      walltime: 48h
      memory: 192GB
      ncpus: 48
      ngpus: 4
    cost_per_hour: 10.0
"""
    
    # Create temporary platform file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(platform_config)
        platform_file = f.name
    
    try:
        # Test with different resource combinations
        from qxub.cli import qxub
        from click.testing import CliRunner
        import logging
        
        # Enable debug logging
        logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
        
        runner = CliRunner()
        
        # Set environment variable for platform path
        platform_dir = os.path.dirname(platform_file)
        os.environ['QXUB_PLATFORM_PATHS'] = platform_dir
        print(f"Platform directory: {platform_dir}")
        print(f"Platform file: {platform_file}")
        
        print("Testing auto-queue selection...")
        
        # Test 1: Small job should select 'small' queue
        print("\n1. Testing small job (2GB, 2 cpus, 30min):")
        result = runner.invoke(qxub, [
            '--queue', 'auto',
            '--resources', 'mem=2GB',
            '--resources', 'ncpus=2', 
            '--resources', 'walltime=30:00',
            '--dry',
            'echo', 'hello'
        ])
        print(f"Exit code: {result.exit_code}")
        if result.output:
            print(f"Output: {result.output}")
        if result.exception:
            print(f"Exception: {result.exception}")
        
        # Test 2: Large memory job should select 'hugemem' queue
        print("\n2. Testing large memory job (500GB):")
        result = runner.invoke(qxub, [
            '--queue', 'auto',
            '--resources', 'mem=500GB',
            '--dry',
            'echo', 'hello'
        ])
        print(f"Exit code: {result.exit_code}")
        if result.output:
            print(f"Output: {result.output}")
            
        # Test 3: GPU job should select 'gpu' queue
        print("\n3. Testing GPU job (2 GPUs):")
        result = runner.invoke(qxub, [
            '--queue', 'auto', 
            '--resources', 'ngpus=2',
            '--dry',
            'echo', 'hello'
        ])
        print(f"Exit code: {result.exit_code}")
        if result.output:
            print(f"Output: {result.output}")
    
    finally:
        # Cleanup
        if os.path.exists(platform_file):
            os.unlink(platform_file)
        # Restore environment
        if 'QXUB_PLATFORM_PATHS' in os.environ:
            del os.environ['QXUB_PLATFORM_PATHS']

if __name__ == "__main__":
    test_auto_queue_selection()