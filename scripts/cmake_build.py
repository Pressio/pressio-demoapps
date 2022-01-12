#!/usr/bin/env python3

import os
import sys
try:
    project_path = f"{os.sep}".join(os.path.abspath(__file__).split(os.sep)[:-2])
    sys.path.append(project_path)
    print(f"Project path was added to sys.path: {project_path}")
except Exception as e:
    print(f"Can not add project path to system path! Exiting!\nERROR: {e}")
    exit(1)

import argparse
import platform
import shutil
import subprocess


class CMakeBuild:
    """ Class responsible fo building library for pressiodemoapps. """
    def __init__(self, build_mode: str = 'Release', enable_omp: bool = False, build_args: list = None):
        # Options
        self.build_mode = build_mode
        self.enable_omp = enable_omp
        # Directories
        self.build_dir = os.path.join(project_path, 'build')
        self._check_and_create_dir(directory=self.build_dir)
        self.temp_dir = os.path.join(self.build_dir, 'temp')
        self._check_and_create_dir(directory=self.temp_dir)
        self.library = os.path.join(self.build_dir, 'pressiodemoapps')
        self._check_and_create_dir(directory=self.library)
        # Build args
        if build_args is not None:
            self.build_args = build_args
        else:
            self.build_args = list()
        self._copy_init_file_to_package()
        self._cli()
        self._print_info()

    def _cli(self):
        """ Support for common line arguments. """
        parser = argparse.ArgumentParser()
        parser.add_argument("--openmp", action="store_true", help="Enables OpenMP if given")
        parser.add_argument("--build_mode", help="Defines build mode: Release or Debug")
        args = parser.parse_args()
        if args.openmp:
            self.enable_omp = True
        if args.openmp == 'Debug':
            self.build_mode = 'Debug'
        else:
            self.build_mode = 'Release'

    @staticmethod
    def _check_and_create_dir(directory: str):
        """ Checks if directory exists. If it does not, creates the path recursively"""
        if not os.path.exists(directory):
            os.makedirs(directory)

    def _copy_init_file_to_package(self):
        """ Copy init file to the library directory """
        file_src = os.path.join(project_path, 'pressiodemoapps', '__init__.py')
        file_dst = os.path.join(self.library, '__init__.py')
        shutil.copyfile(file_src, file_dst)

    def _print_info(self):
        """ Prints information regarding paths and options. """
        print('\n')
        print('======> PROJECT DIRECTORIES <======')
        print(f'===> Project path: {project_path}')
        print(f'===> Build directory: {self.build_dir}')
        print(f'===> Temp directory: {self.temp_dir}')
        print(f'===> Library directory: {self.library}')
        print('\n')
        print('======> PROJECT OPTIONS <======')
        print(f'===> Project build mode: {self.build_mode}')
        print(f'===> Enable OpenMP: {self.enable_omp}')
        print('\n')

    def _check_os_env_vars(self):
        """ Checks environment variables. """
        if "CXX" not in os.environ:
            if platform.system() != 'Windows':
                compiler = subprocess.run(['which', 'g++'], capture_output=True)
                if compiler.returncode == 1:
                    print(f"CXX env var missing, needs to point to your target C++ compiler")
                    exit(1)
                else:
                    cxx = compiler.stdout.decode('utf-8').replace('\n', '')
                    os.environ["CXX"] = cxx
            else:
                print(f"CXX env var missing, needs to point to your target C++ compiler")
                exit(1)
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            self.build_args.append('-j4')

    def build(self):
        """ Builds library with Cmake. """
        self._check_os_env_vars()
        cmake_args = [f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={self.library}", f"-DPYTHON_EXECUTABLE={sys.executable}",
                      f"-DCMAKE_BUILD_TYPE={self.build_mode}", "-DPRESSIODEMOAPPS_ENABLE_BINDINGS=On",
                      "-DCMAKE_VERBOSE_MAKEFILE=On", f"-DPRESSIODEMOAPPS_ENABLE_OPENMP={self.enable_omp}"]

        subprocess.check_call(["cmake", project_path] + cmake_args, cwd=self.temp_dir)
        subprocess.check_call(["cmake", "--build", "."] + self.build_args, cwd=self.temp_dir)


if __name__ == '__main__':
    CMakeBuild().build()
