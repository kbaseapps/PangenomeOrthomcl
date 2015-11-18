FROM kbase/sdkbase:latest
MAINTAINER KBase Developer
# -----------------------------------------

# Insert apt-get instructions here to install
# any required dependencies for your module.

# RUN apt-get update
RUN apt-get update && \
  echo 'mysql-server mysql-server/root_password password 12345' | debconf-set-selections && \
  echo 'mysql-server mysql-server/root_password_again password 12345' | debconf-set-selections && \
  apt-get -y install mysql-server
RUN service mysql start
WORKDIR /tmp
RUN wget http://www.micans.org/mcl/src/mcl-12-068.tar.gz
RUN tar xzf mcl-12-068.tar.gz
RUN wget http://orthomcl.org/common/downloads/software/v2.0/orthomclSoftware-v2.0.9.tar.gz
RUN tar xzf orthomclSoftware-v2.0.9.tar.gz
WORKDIR /tmp/mcl-12-068
RUN ./configure --prefix=/usr
RUN make install
WORKDIR /tmp/orthomclSoftware-v2.0.9
RUN cp ./bin/* /kb/deployment/plbin/
RUN cp -r ./lib/perl/* /kb/deployment/lib/

RUN apt-get -y install mc
# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work

WORKDIR /kb/module

RUN make

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
