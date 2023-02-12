from setuptools import setup, find_packages

setup(name='tetres',
      version='1.0.0',
      description='A package for discrete time trees',
      url='https://github.com/bioDS/Centroid-Code',  # this is only for the paper release
      author='Lars Berling',
      author_email='berlinglars96@gmail.com',
      license='MIT',
      packages=find_packages(),
      install_requires=['ete3',
                        'pandas',
                        'seaborn',
                        'rpy2'
                        ],
      zip_safe=False)