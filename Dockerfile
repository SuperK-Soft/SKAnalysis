### Created by Dr. Benjamin Richards

### Download base image from repo
FROM toolframework/centos7

### Run the following commands as super user (root):
USER root

Run git clone https://github.com/ToolFramework/ToolApplication.git

### Open terminal
CMD ["/bin/bash"]