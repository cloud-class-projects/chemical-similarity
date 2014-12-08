sudo mkdir -p /data/shardb
sudo chown ubuntu /data
sudo chown ubuntu /data/shardb
sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 7F0CEB10
echo 'deb http://downloads-distro.mongodb.org/repo/ubuntu-upstart dist 10gen' | sudo tee /etc/apt/sources.list.d/mongodb.list
sudo apt-get update
sudo apt-get install -y mongodb-org
ssh-keygen -b 4096 -t rsa -f /home/ubuntu/.ssh/id_rsa -P ""