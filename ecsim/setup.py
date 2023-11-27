from setuptools import setup,find_packages

setup(name='ecsim',
      version = '1.0.0',
      description = 'Circular DNA simualation tools',
      author='LiuKe',
      url='https://github.com/qukunLab/ecsim',
      packages=find_packages(),
      install_requires=[],
      package_dir={'ecsim': 'ecsim'},
      package_data={'ecsim/resources': ['ecsim/resources/pbsim2/*model'], 'ecsim/resources': ['ecsim/resources/profile/*bed']},
      entry_points={
          	'console_scripts': ['ecsim = ecsim.ecsim:main' ,'downsample = ecsim.downsample:main'],
      },
      classifiers=[
          'License :: OSI Approved :: MIT License'
      ],
     )

