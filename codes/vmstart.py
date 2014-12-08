###############################################################################################
###  Creating and starting the VMs  
###  Abhik Seal
###  1 Nov 2014
###############################################################################################
import cloudmesh
import sys
import os

# Creating cloud on India server
def createVM(N):
    username = cloudmesh.load().username()
    print username
    mesh = cloudmesh.mesh("mongo")
    mesh.activate(username)
    cloudmesh.shell("cloud on india")
    flavor = mesh.flavor('india', 'm1.large')
    image=mesh.image('india','futuregrid/ubuntu-14.04')
    #vm_ip=[]
    for i in range(0,N):
        result = mesh.start(cloud='india',
                        cm_user_id=username,
                        flavor=flavor,
                        image=image)
        server = result['server']['id']
        ip=mesh.assign_public_ip('india', server, username)
        print ip
        #vm_ip.append(ip)
        try:
            result = mesh.wait(ipaddr=ip, command="ls -al", interval=10, retry=5)
            cloudmesh.banner("INSTALLING MONGO TO THE VM  "+str(ip))
            result=mesh.ssh_execute(username='ubuntu',ipaddr=ip,command='ls -al',pkey=None)
            print result
            script ="scp  filescript.sh  ubuntu@"+ip+":/home/ubuntu/"
            print script
            os.system(script)
            result=mesh.ssh_execute(username='ubuntu',ipaddr=ip,command='bash -s < ./filescript.sh',pkey=None)
            print result
            result=mesh.ssh_execute(username='ubuntu',ipaddr=ip,command='ls -al',pkey=None)
            print result
        except:
            print "Authentication failed when connecting to %s" % ip
            sys.exit(1)
    print "All VMs Mongo Installed.. "
    #return(vm_ip)



if __name__ == '__main__':
    
    createVM(N=3)
