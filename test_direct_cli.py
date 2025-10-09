#!/usr/bin/env python3
"""
Direct test of CLI auto-selection with manual monkeypatching
"""

import os
import sys

# Add the parent directory to Python path  
sys.path.insert(0, '/g/data/a56/software/qsub_tools')

def test_direct_cli_auto_selection():
    """Test by directly calling the CLI auto-selection logic"""
    
    # Mock the parameters that would come from Click
    params = {
        "queue": "auto",
        "resources": ["mem=2GB", "ncpus=2", "walltime=30:00"],
        "name": "test-job",
        "project": "a56",
        "joblog": None,
        "out": "/scratch/a56/jr9959/qt/test/out",
        "err": "/scratch/a56/jr9959/qt/test/err",
        "dry": True,
        "quiet": False
    }
    
    # Set the environment variable
    os.environ['QXUB_PLATFORM_PATHS'] = '/g/data/a56/software/qsub_tools/test_platforms'
    
    # Test the auto-selection code directly
    if params["queue"] == "auto":
        try:
            from qxub.platform import PlatformLoader
            from qxub.resource_utils import parse_memory_size, parse_walltime
            from pathlib import Path
            
            # Check for QXUB_PLATFORM_PATHS environment variable
            platform_paths_env = os.environ.get('QXUB_PLATFORM_PATHS')
            print(f"Environment variable QXUB_PLATFORM_PATHS: {platform_paths_env}")
            
            if platform_paths_env:
                search_paths = [Path(p.strip()) for p in platform_paths_env.split(':')]
                loader = PlatformLoader(search_paths=search_paths)
                print(f"Created loader with search paths: {search_paths}")
            else:
                loader = PlatformLoader()
                print("Created loader with default search paths")
            
            platform_names = loader.list_platforms()
            print(f"Available platforms: {platform_names}")
            
            if not platform_names:
                print("No platforms available for auto queue selection, using 'normal'")
                params["queue"] = "normal"
            else:
                # Build requirements from resources
                requirements = {}
                if params.get("resources"):
                    for resource in params["resources"]:
                        if "=" in resource:
                            key, value = resource.split("=", 1)
                            if key == "mem":
                                requirements["memory"] = value  # Keep as string for condition parsing
                            elif key == "walltime":
                                requirements["walltime"] = value  # Keep as string for condition parsing
                            elif key == "ncpus":
                                requirements["cpus"] = int(value)  # Use "cpus" not "ncpus"
                            elif key in ["ngpus", "gpu"]:
                                gpu_count = int(value) if value.isdigit() else 1
                                requirements["gpu_requested"] = gpu_count  # Use "gpu_requested" for conditions
                                requirements["gpus"] = gpu_count  # Keep "gpus" for validation
                
                print(f"Built requirements: {requirements}")
                
                # Try to find best queue from any platform
                best_queue = None
                best_cost = float('inf')
                
                for platform_name in platform_names:
                    platform = loader.get_platform(platform_name)
                    if not platform:
                        continue
                        
                    try:
                        selected_queue = platform.select_queue(requirements)
                        print(f"Platform {platform_name} selected queue: {selected_queue}")
                        if selected_queue:
                            # Get estimated cost for comparison
                            queue = platform.get_queue(selected_queue)
                            if queue:
                                # Estimate cost based on cores and walltime
                                cores = requirements.get("cpus", 1)
                                walltime_hours = requirements.get("walltime", 3600)
                                # Convert walltime to hours if it's a string
                                if isinstance(walltime_hours, str):
                                    from qxub.resource_utils import parse_walltime
                                    walltime_hours = parse_walltime(walltime_hours) or 1.0
                                elif isinstance(walltime_hours, (int, float)):
                                    walltime_hours = walltime_hours / 3600.0  # Convert seconds to hours
                                
                                cost = queue.estimate_su_cost(cores, walltime_hours)
                                print(f"Cost for queue {selected_queue}: {cost}")
                                if cost < best_cost:
                                    best_cost = cost
                                    best_queue = selected_queue
                    except Exception as e:
                        print(f"Failed to select queue from platform {platform.name}: {e}")
                        continue
                
                if best_queue:
                    params["queue"] = best_queue
                    print(f"Auto-selected queue: {best_queue}")
                else:
                    print("No suitable queue found for requirements, using 'normal'")
                    params["queue"] = "normal"
                    
        except ImportError:
            print("Platform system not available for auto queue selection, using 'normal'")
            params["queue"] = "normal"
        except Exception as e:
            print(f"Auto queue selection failed: {e}, using 'normal'")
            params["queue"] = "normal"
    
    print(f"Final queue selection: {params['queue']}")

if __name__ == "__main__":
    test_direct_cli_auto_selection()