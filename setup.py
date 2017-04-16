from setuptools import setup
from setuptools import find_packages

          
setup(name='spectrumfitter',
      version='0.1.0',
      description='basic utility for fitting dielectric spectra',
      author='Daniel C Elton',
      author_email='delton17@gmail.com',
      download_url='https://github.com/delton137/spectrumfitter/archive/master.zip',
      url='https://github.com/delton137/spectrumfitter',
      license='MIT',
      install_requires=['numpy', 'matplotlib', 'scipy'],
      packages=find_packages(),
      zip_safe=False)

