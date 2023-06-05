
## Instructions for running Autogluon (Python scripts)


### pull docker image autogluon from registry
```
docker pull autogluon/autogluon:0.6.1-cpu-jupyter-ubuntu20.04-py3.8
```

* Digest: sha256:ace9e4217b43d9072ca9087d3fed819a5c4a0c01ec9044e8fde4a0692f0e1f6e
* Status: Downloaded newer image for autogluon/autogluon:0.6.1-cpu-jupyter-ubuntu20.04-py3.8
* docker.io/autogluon/autogluon:0.6.1-cpu-jupyter-ubuntu20.04-py3.8

### start container interactively (in screen) and run script
```
docker run --rm -it -v $(pwd):/result autogluon/autogluon:0.6.1-cpu-jupyter-ubuntu20.04-py3.8 /bin/bash
cd result
python3 <script>.py
```

