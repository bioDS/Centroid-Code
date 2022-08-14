from setuptools import setup, find_packages

setup(name='treeoclock',
      version='1.0.0',
      description='A package for discrete time trees',
      url='https://github.com/bioDS/Centroid-Code',  # this is only for the paper release
      author='Lars Berling',
      author_email='berlinglars96@gmail.com',
      license='MIT',
      packages=find_packages(),
      install_requires=['ete3',
                        'pandas',
                        'seaborn'
                        ],
      zip_safe=False)