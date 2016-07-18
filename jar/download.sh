#!/usr/bin/env bash

JAR_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cd ${JAR_DIR}

echo "==> jrunlist"
wget http://repo1.maven.org/maven2/com/github/egateam/jrunlist/0.1.4/jrunlist-0.1.4-jar-with-dependencies.jar \
    -O jrunlist.jar

echo "==> jrange"
wget -N http://repo1.maven.org/maven2/com/github/egateam/jrange/0.0.3/jrange-0.0.3-jar-with-dependencies.jar \
    -O jrange.jar
