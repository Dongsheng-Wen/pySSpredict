from setuptools import setup

#name1 = 'sspredict_local.plot = sspredict.plot_ss_ternary:main'
name1 = 'sspredict.predict = sspredict.master:main'

setup(
    name='SSPredict',
    version='v1.1.0',
    license ='MIT',
    author='Dongsheng Wen, Michael S. Titus',
    author_email='wen94@purdue.edu, titus9@purdue.edu',
    description='Python-based Solid-Solution Strengthening Prediction Tool',
    packages = ['sspredict','sspredict/make_prediction'],
    entry_points={
          'console_scripts': [name1]
      },
    platforms='any',
    install_requires=[
    'numpy >= 1.13.3',
    'matplotlib >= 2.1.0',
    'pandas >= 0.23.4',
    'python-ternary >= 1.0.8',
    'scipy >= 1.5.1']
)
