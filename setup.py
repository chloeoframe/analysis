from setuptools import setup, find_packages
import sys

try:
    import mdtraj
except ImportError:
    print('Building and running analysis requires mdtraj. See '
          'http://mdtraj.org/latest/installation.html for help!')
    sys.exit(1)

setup(name='analysis',
      version='0.1',
      description=('Analysis scripts for molecular dynamics simulations' +
                  ' of stratum corneum lipid multilayers.'),
      url='https://github.com/uppittu11/analysis',
      author='Parashara Shamaprasad',
      author_email='p.shama@vanderbilt.edu',
      license='MIT',
      packages=find_packages(),
      package_dir={'analysis': 'analysis'},
      include_package_data=True,
      install_requires=["mdtraj"],
      entry_points={
          "console_scripts" : ["analyze=bin.analyze:main"],
          },
      scripts=['bin/progress.sh'],
)
