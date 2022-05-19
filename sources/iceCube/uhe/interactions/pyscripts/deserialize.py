import javaobj
import numpy as np
import click
from tqdm import tqdm
import sys

@click.command()
@click.option("--javafile",required=True)
@click.option("--material",type=int,required=True)
def deserialize(javafile,material):
    intMtx = getIntMtx(javafile)
    print(intMtx)
    transMtx  = intMtx.transferMtx
    survivMtx = intMtx.transferAMtx
    sigmaMtx  = intMtx.sigmaMtx
    ydsdyMtx  = intMtx.inelasticityMtx

    path = ['iceCube/uhe/interactions/ice/new/','iceCube/uhe/interactions/rock/new/']
    print(transMtx[699])
    filenamestr = javafile.split('/')
    filename = filenamestr[-1]
    with open(f'{path[material]}{filename}_sigma.txt','w') as f:
        for i,sigma in enumerate(sigmaMtx):
            f.write(str(sigma))
            if i+1 < len(sigmaMtx):
                f.write(',')
    with open(f'{path[material]}{filename}_inela.txt','w') as f:
        for i,inela in enumerate(ydsdyMtx):
            f.write(str(inela))
            if i+1 < len(ydsdyMtx):
                f.write(',')
    with open(f'{path[material]}{filename}_trans.txt','w') as f:
        for i,transs in enumerate(transMtx):
            for j,trans in enumerate(transs):
                f.write(str(trans))
                if j+1 < len(transs):
                    f.write(',')
                elif i+1 < len(transMtx):
                    f.write('\n')
                else:
                    pass
    with open(f'{path[material]}{filename}_surviv.txt','w') as f:
        for i,survivs in enumerate(survivMtx):
            for j,surviv in enumerate(survivs):
                f.write(str(surviv))
                if j+1 < len(survivs):
                    f.write(',')
                elif i+1 < len(survivMtx):
                    f.write('\n')
                else:
                    pass

def getIntMtx(javafile):
    with open(javafile,"rb") as f:
        sobj = f.read()

    intMtx = javaobj.loads(sobj)
    return intMtx

if __name__ == "__main__":
    deserialize()
