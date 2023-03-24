from setuptools import setup, find_packages

setup(name='tetres',
      version='0.0.1',
      description='TimE TREe Statistics package',
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