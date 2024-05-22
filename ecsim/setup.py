from setuptools import setup,find_packages

setup(
    name='ecsim',
    version='1.0.1',
    description='Circular DNA simulation tools',
    author='Ke Liu',
    url='https://github.com/qukunLab/eccDNABenchmarking',
    packages=find_packages(),
    package_dir={'ecsim': 'ecsim'},
    install_requires=[
    ],
    include_package_data=True,
    package_data={
        'ecsim': ['resource/pbsim2/*model'],
    },
    entry_points={
        'console_scripts': [
            'ecsim=ecsim.ecsim:main',
            'downsample=ecsim.downsample:main'
        ],
    },
    classifiers=[
        'License :: OSI Approved :: MIT License',
    ],
)
