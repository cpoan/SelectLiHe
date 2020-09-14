#!/usr/bin/env python

import sys,os
if __name__ == "__main__":
    dataPeriod = "P17B"
    RunType = ""
    program = "main"

    site = "EH1"
    GoodRunList = "/dybfs/users/chenpoan/"+dataPeriod+"/"+site+"/"
    AppDir = "/dybfs/users/chenpoan/EventSelect/SelectLiHe/"
    OutPutDir = AppDir+"rootFile/"
    scriptDir = AppDir+"tcshFile/"

    ListOfList = os.listdir( GoodRunList )
    for run in ListOfList:
        listName = run.split('.')[0]
        runNo = run.split('.')[0]
        print runNo
        base = "p_"+runNo
        cshfile = scriptDir+base+".tcsh"
        FILE = open(cshfile,"w")
        FILE.writelines("#!/bin/tcsh \n")
        FILE.writelines("rerun:\n")
        FILE.writelines("cd "+AppDir+" \n")
        FILE.writelines("source  /afs/ihep.ac.cn/soft/dayabay/NuWa-slc6/opt/external/ROOT/5.26.00e_python2.7/x86_64-slc6-gcc44-opt/root/bin/thisroot.sh\n")
        FILE.writelines("./"+program+"1 "+GoodRunList+listName+".list "+OutPutDir+base+".root \n")
        FILE.writelines("if(\"$?\" == \"0\")then\n")
        FILE.writelines("\techo \"Let's do it\"\n")
        FILE.writelines("else\n")
        FILE.writelines("\techo \"Resource again\"\n")
        FILE.writelines("\tgoto rerun\n")
        FILE.writelines("fi")
        FILE.close()
        os.system("chmod a+x "+cshfile)
        os.system( "/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin/hep_sub "+cshfile+" -g dyw" )

    site = "EH2"
    GoodRunList = "/dybfs/users/chenpoan/"+dataPeriod+"/"+site+"/"
    ListOfList = os.listdir( GoodRunList )
    for run in ListOfList:
        listName = run.split('.')[0]
        runNo = run.split('.')[0]
        print runNo
        base = "p_"+runNo
        cshfile = scriptDir+base+".tcsh"
        FILE = open(cshfile,"w")
        FILE.writelines("#!/bin/tcsh \n")
        FILE.writelines("rerun:\n")
        FILE.writelines("cd "+AppDir+" \n")
        FILE.writelines("source  /afs/ihep.ac.cn/soft/dayabay/NuWa-slc6/opt/external/ROOT/5.26.00e_python2.7/x86_64-slc6-gcc44-opt/root/bin/thisroot.sh\n")
        FILE.writelines("./"+program+"2 "+GoodRunList+listName+".list "+OutPutDir+base+".root \n")
        FILE.writelines("if(\"$?\" == \"0\")then\n")
        FILE.writelines("\techo \"Let's do it\"\n")
        FILE.writelines("else\n")
        FILE.writelines("\techo \"Resource again\"\n")
        FILE.writelines("\tgoto rerun\n")
        FILE.writelines("fi")
        FILE.close()
        os.system("chmod a+x "+cshfile)
        os.system( "/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin/hep_sub "+cshfile+" -g dyw" )

    site = "EH3"
    GoodRunList = "/dybfs/users/chenpoan/"+dataPeriod+"/"+site+"/"
    ListOfList = os.listdir( GoodRunList )
    for run in ListOfList:
        listName = run.split('.')[0]
        runNo = run.split('.')[0]
        print runNo
        base = "p_"+runNo
        cshfile = scriptDir+base+".tcsh"
        FILE = open(cshfile,"w")
        FILE.writelines("#!/bin/tcsh \n")
        FILE.writelines("rerun:\n")
        FILE.writelines("cd "+AppDir+" \n")
        FILE.writelines("source  /afs/ihep.ac.cn/soft/dayabay/NuWa-slc6/opt/external/ROOT/5.26.00e_python2.7/x86_64-slc6-gcc44-opt/root/bin/thisroot.sh\n")
        FILE.writelines("./"+program+"3 "+GoodRunList+listName+".list "+OutPutDir+base+".root \n")
        FILE.writelines("if(\"$?\" == \"0\")then\n")
        FILE.writelines("\techo \"Let's do it\"\n")
        FILE.writelines("else\n")
        FILE.writelines("\techo \"Resource again\"\n")
        FILE.writelines("\tgoto rerun\n")
        FILE.writelines("fi")
        FILE.close()
        os.system("chmod a+x "+cshfile)
        os.system( "/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin/hep_sub "+cshfile+" -g dyw" )
