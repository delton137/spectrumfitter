from setuptools import setup
from setuptools import find_packages

__version__ = '0.1.0'

def readme():
    with open('README.md') as f:
        return f.read()
          
setup(name='spectrumfitter',
      version='0.1.0',
      description='basic utility for fitting dielectric spectra',
      author='Daniel C Elton',
      author_email='delton17@gmail.com',
      download_url='https://github.com/delton137/spectrumfitter/archive/master.zip',
      url='https://github.com/delton137/spectrumfitter',
      include_package_data=True,
      classifiers=[
        'Topic :: Scientific/Engineering :: Physics',
	'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
      	'Intended Audience :: Science/Research'
	],
      license='MIT',
      install_requires=['numpy', 'matplotlib', 'scipy'],
      packages=find_packages(),
      zip_safe=False)

