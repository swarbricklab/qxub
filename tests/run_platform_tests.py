#!/usr/bin/env python3
"""
Comprehensive test runner for qxub platform system.

Runs all test modules and provides a summary report.
"""

import sys
import subprocess
import time
from pathlib import Path

def run_test_script(script_path: Path) -> tuple[int, str, float]:
    """
    Run a test script and return exit code, output, and execution time.
    
    Returns:
        Tuple of (exit_code, output, execution_time)
    """
    start_time = time.time()
    
    try:
        result = subprocess.run(
            [sys.executable, str(script_path)],
            capture_output=True,
            text=True,
            timeout=60  # 60 second timeout
        )
        execution_time = time.time() - start_time
        return result.returncode, result.stdout + result.stderr, execution_time
    
    except subprocess.TimeoutExpired:
        execution_time = time.time() - start_time
        return -1, "Test timed out after 60 seconds", execution_time
    
    except Exception as e:
        execution_time = time.time() - start_time
        return -2, f"Failed to run test: {str(e)}", execution_time

def main():
    """Run all test scripts and provide summary."""
    print("🚀 Running qxub platform system test suite...\n")
    
    test_dir = Path(__file__).parent
    test_scripts = [
        ("Resource Utilities", test_dir / "test_resource_utils.py"),
        ("Platform System", test_dir / "test_platform.py"),
        ("Configuration", test_dir / "test_config.py"),
    ]
    
    results = []
    total_time = 0
    
    for test_name, script_path in test_scripts:
        print(f"{'='*60}")
        print(f"Running {test_name} tests...")
        print(f"{'='*60}")
        
        if not script_path.exists():
            print(f"❌ Test script not found: {script_path}")
            results.append((test_name, -3, "Script not found", 0))
            continue
        
        exit_code, output, exec_time = run_test_script(script_path)
        total_time += exec_time
        
        # Print the test output
        print(output)
        
        # Store result
        if exit_code == 0:
            status = "✅ PASSED"
        elif exit_code == -1:
            status = "⏰ TIMEOUT"
        elif exit_code == -2:
            status = "💥 ERROR"
        else:
            status = "❌ FAILED"
        
        results.append((test_name, exit_code, status, exec_time))
        print(f"{status} ({exec_time:.2f}s)\n")
    
    # Print summary
    print(f"{'='*60}")
    print("📊 TEST SUMMARY")
    print(f"{'='*60}")
    
    passed = 0
    failed = 0
    
    for test_name, exit_code, status, exec_time in results:
        print(f"{test_name:<25} {status} ({exec_time:.2f}s)")
        if exit_code == 0:
            passed += 1
        else:
            failed += 1
    
    print(f"\n📈 Results:")
    print(f"   Passed: {passed}")
    print(f"   Failed: {failed}")
    print(f"   Total:  {len(results)}")
    print(f"   Time:   {total_time:.2f}s")
    
    if failed == 0:
        print("\n🎉 All tests passed! Platform system is ready for integration.")
        return 0
    else:
        print(f"\n❌ {failed} test group(s) failed. Please review the output above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())