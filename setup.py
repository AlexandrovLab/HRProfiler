from setuptools import setup, find_packages
import os
import shutil


#remove the dist folder first if exists
if os.path.exists("dist"):
  shutil.rmtree("dist")

def readme():
  this_directory = os.path.abspath(os.path.dirname(__file__))
  with open(os.path.join(this_directory, 'README.md'), encoding='latin-1') as f:
    long_description = f.read()
    return(long_description)
    
VERSION = '0.0.2'

def write_version_py(filename='HRProfiler/version.py'):
  # Copied from numpy setup.py
  cnt = """
# THIS FILE IS GENERATED FROM HRProfiler SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
  
  """
  fh = open(filename, 'w')
  fh.write(cnt % {'version': VERSION,})
  fh.close()

write_version_py()

setup(
    name="HRProfiler",
    version=VERSION,
    description="HRD Prediction from WGS and WES samples",
    long_description=readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/AlexandrovLab/HRProfiler.git",
    author="Ammal Abbasi",
    author_email="amabbasi@ucsd.edu",
    license="UCSD",
    packages=["HRProfiler", "HRProfiler.scripts"],
    python_requires=">=3.9",
    install_requires=[
        "sigProfilerPlotting>=1.4.1",
        "SigProfilerMatrixGenerator>=1.3.5",
        "joblib>=0.16.0",
        "scikit-plot==0.3.7",
        "scikit-learn>=1.1.3",
         "seaborn",
         "numpy >= 2.0.0",
         "pandas >= 2.0.0",
         "matplotlib >= 3.4.3",

    ],
    include_package_data=True,
    zip_safe=False,
)

