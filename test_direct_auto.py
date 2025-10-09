#!/usr/bin/env python3
"""
Direct test of auto-queue selection logic
"""

import tempfile
import os
import sys

# Add the parent directory to Python path
sys.path.insert(0, '/g/data/a56/software/qsub_tools')

def test_direct_auto_selection():
    """Test the auto-queue selection logic directly"""
    
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
    import tempfile
    with tempfile.TemporaryDirectory() as temp_dir:
        platform_file = os.path.join(temp_dir, "test_platform.yaml")
        with open(platform_file, 'w') as f:
            f.write(platform_config)
    
        try:
            # Test the platform loading and queue selection directly
            from qxub.platform import PlatformLoader
            from qxub.resource_utils import parse_memory_size, parse_walltime
            import logging
            
            # Enable debug logging
            logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
            
            # Set platform path
            platform_dir = temp_dir
            os.environ['QXUB_PLATFORM_PATHS'] = platform_dir
            
            print(f"Platform file: {platform_file}")
            print(f"Platform directory: {platform_dir}")
            print(f"Environment QXUB_PLATFORM_PATHS: {os.environ.get('QXUB_PLATFORM_PATHS')}")
            
            # Verify file exists
            print(f"File exists: {os.path.exists(platform_file)}")
            
            # Load platforms
            from pathlib import Path
            loader = PlatformLoader(search_paths=[Path(platform_dir)])
            platform_names = loader.list_platforms()
            print(f"Available platforms: {platform_names}")
            
            if platform_names:
                platform = loader.get_platform(platform_names[0])
                print(f"Platform: {platform.name}")
                print(f"Queues: {platform.list_queues()}")
                
                # Test 1: Small job requirements
                print("\n1. Testing small job requirements:")
                requirements = {
                    "memory": parse_memory_size("2GB"),
                    "ncpus": 2,
                    "walltime": parse_walltime("30:00")
                }
                print(f"Requirements: {requirements}")
                
                selected_queue = platform.select_queue(requirements)
                if selected_queue:
                    print(f"Selected queue: {selected_queue}")
                else:
                    print("No queue selected")
                
                # Test 2: Large memory requirements
                print("\n2. Testing large memory requirements:")
                requirements = {
                    "memory": parse_memory_size("500GB")
                }
                print(f"Requirements: {requirements}")
                
                selected_queue = platform.select_queue(requirements)
                if selected_queue:
                    print(f"Selected queue: {selected_queue}")
                else:
                    print("No queue selected")
                
                # Test 3: GPU requirements
                print("\n3. Testing GPU requirements:")
                requirements = {
                    "ngpus": 2
                }
                print(f"Requirements: {requirements}")
                
                selected_queue = platform.select_queue(requirements)
                if selected_queue:
                    print(f"Selected queue: {selected_queue}")
                else:
                    print("No queue selected")
            else:
                print("No platforms loaded!")
        
        except Exception as e:
            print(f"Error: {e}")
            import traceback
            traceback.print_exc()

if __name__ == "__main__":
    test_direct_auto_selection()