#!/usr/bin/env python3
"""
Simple test of auto-queue selection directly in CLI context
"""

import os
import sys
import tempfile

# Add the parent directory to Python path
sys.path.insert(0, '/g/data/a56/software/qsub_tools')

def test_simple_auto_selection():
    """Test auto-selection logic by calling qxub directly with our test platform"""
    from qxub.cli import qxub
    from click.testing import CliRunner
    
    runner = CliRunner()
    
    # Use our persistent test platform
    platform_dir = '/g/data/a56/software/qsub_tools/test_platforms'
    
    print("Testing auto-queue selection in isolated CLI environment...")
    
    # Test 1: Small job with dry run
    print("\n1. Testing small job with auto queue selection:")
    result = runner.invoke(qxub, [
        '--queue', 'auto',
        '--resources', 'mem=2GB',
        '--resources', 'ncpus=2', 
        '--resources', 'walltime=30:00',
        '--dry',
        'echo', 'hello'
    ], env={'QXUB_PLATFORM_PATHS': platform_dir})
    
    print(f"Exit code: {result.exit_code}")
    if result.output:
        print(f"Output: {result.output}")
    if result.exception:
        print(f"Exception: {result.exception}")
        import traceback
        traceback.print_exception(type(result.exception), result.exception, result.exception.__traceback__)

    # Test 2: Large memory job 
    print("\n2. Testing large memory job with auto queue selection:")
    result = runner.invoke(qxub, [
        '--queue', 'auto',
        '--resources', 'mem=500GB',
        '--dry',
        'echo', 'hello'
    ], env={'QXUB_PLATFORM_PATHS': platform_dir})
    
    print(f"Exit code: {result.exit_code}")
    if result.output:
        print(f"Output: {result.output}")

    # Test 3: GPU job
    print("\n3. Testing GPU job with auto queue selection:")
    result = runner.invoke(qxub, [
        '--queue', 'auto',
        '--resources', 'ngpus=2',
        '--dry',
        'echo', 'hello'
    ], env={'QXUB_PLATFORM_PATHS': platform_dir})
    
    print(f"Exit code: {result.exit_code}")
    if result.output:
        print(f"Output: {result.output}")

if __name__ == "__main__":
    test_simple_auto_selection()