#!/bin/bash
# bash_master.sh: ssh'ing into corn servers and executing MA code

echo "Please enter password for lanemc@corn.stanford.edu"
read PASSWORD


#cd ~/Git/masters-thesis/run_these/
for i in 'ls ~/Git/masters-thesis/run_these/';
# ssh into corn
do sshpass -p $PASSWORD ssh lanemc@corn.stanford.edu 'cd ~/masters/python/; python "$i" &' &;

done


# prepare corn for running a long program in a screen session
#pagsh
#kinit;aklog 
#$PASSWORD
 
#screen -S masters

#keeptoken
#source /tmp/.krbhold_lanemc.csh


#python mastersRun_chunking10e4_brownian_i.py


# to detach screen execute keyboard shortcut
#xdotool key "ctrl+a+d" 
