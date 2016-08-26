from setuptools import setup

setup(name='eqstats',
      version='0.1',
      description='Tools useful for statistical analysis of earthquake catalogs',
      url='http://bitbucket.org/egdaub/eqstats',
      author='Eric Daub',
      author_email='egdaub@memphis.edu',
      packages=['eqstats'],
      install_requires=['numpy', 'scipy'],
      test_suite='nose.collector',
      tests_require=['nose'])
