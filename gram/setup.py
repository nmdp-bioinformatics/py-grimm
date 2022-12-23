from setuptools import setup
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
   name='gramm',
   version='0.0.1',
   description='Imputation haplotypes in nuclear families',
   long_description=long_description,
   long_description_content_type='text/markdown',
   author='YOLO lab',
   author_email='louzouy@math.biu.ac.il',
   packages=['gramm'],
   install_requires=['matplotlib==3.3.2', 'networkx==2.2', 'pandas==1.4.1', 'numpy==1.19.2', 'typing==3.6.4', 'py-ard==0.6.9'],
)