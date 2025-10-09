#!/usr/bin/env python3
"""
Debug auto-queue selection to understand why it always chooses 'normal'
"""

import os
import sys

# Add the parent directory to Python path
sys.path.insert(0, '/g/data/a56/software/qsub_tools')

def debug_auto_selection():
    """Debug the auto-selection logic step by step"""
    
    from qxub.platform import PlatformLoader
    from qxub.resource_utils import parse_memory_size, parse_walltime
    from pathlib import Path
    import logging
    
    # Enable debug logging
    logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
    
    # Set up test platform
    platform_dir = '/g/data/a56/software/qsub_tools/test_platforms'
    
    print(f"Loading platforms from: {platform_dir}")
    
    # Load platforms directly with search paths
    loader = PlatformLoader(search_paths=[Path(platform_dir)])
    platform_names = loader.list_platforms()
    print(f"Available platforms: {platform_names}")
    
    if not platform_names:
        print("No platforms found!")
        return
        
    platform = loader.get_platform(platform_names[0])
    print(f"Platform: {platform.name}")
    print(f"Available queues: {platform.list_queues()}")
    
    # Test different requirements sets
    test_cases = [
        {
            "name": "Small job",
            "requirements": {
                "memory": "2GB",  # Keep as string for validation
                "cpus": 2,
                "walltime": "30:00"
            }
        },
        {
            "name": "Large memory job", 
            "requirements": {
                "memory": "500GB"  # Keep as string for validation
            }
        },
        {
            "name": "GPU job",
            "requirements": {
                "gpus": 2  # Note: using "gpus" not "ngpus"
            }
        }
    ]
    
    # Test different requirements sets
    test_cases = [
        {
            "name": "Small job",
            "requirements": {
                "memory": "2GB",
                "cpus": 2,
                "walltime": "30:00"
            }
        },
        {
            "name": "Large memory job", 
            "requirements": {
                "memory": "500GB"
            }
        },
        {
            "name": "GPU job",
            "requirements": {
                "gpus": 2  # For validation, but need gpu_requested for condition
            }
        }
    ]
    
    for test_case in test_cases:
        print(f"\n--- Testing: {test_case['name']} ---")
        print(f"Requirements: {test_case['requirements']}")
        
        # For GPU jobs, we need to add gpu_requested for the condition parser
        selection_requirements = test_case['requirements'].copy()
        if 'gpus' in selection_requirements:
            selection_requirements['gpu_requested'] = selection_requirements['gpus']
        
        # Test queue selection
        selected_queue = platform.select_queue(selection_requirements)
        print(f"Selected queue: {selected_queue}")
        
        # Check each queue individually (use original requirements for validation)
        for queue_name in platform.list_queues():
            queue = platform.get_queue(queue_name)
            print(f"  Queue {queue_name}: limits={queue.limits}")
            
            # Check if this queue can handle the requirements
            validation_result = queue.validate_resources(test_case['requirements'])
            print(f"    Validation: {validation_result}")

if __name__ == "__main__":
    debug_auto_selection()