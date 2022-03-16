from setuptools import setup, find_packages

with open('README.md', 'r') as readme_file:
    readme = readme_file.read()

requirements = ['rdkit-pypi>=2021']

setup(
    name='magicmethyl',
    version='0.0.1',
    author='Sam Mun',
    author_email='therealsam@berkeley.edu',
    description='A package to find conformationally unique "magic methyls"',
    long_description=readme,
    long_description_content_type='text/markdown',
    url='https://github.com/moleculeclub/magicmethyl',
    packages=find_packages(),
    install_requires=requirements,
    classifiers=[
        'Programming Language :: Python :: 3.9',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)