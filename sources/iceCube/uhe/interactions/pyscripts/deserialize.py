import javaobj
import numpy as np
import click

@click.command()
@click.option("--javafile",required=True)
def deserialize(javafile):
    intMtx = getIntMtx(javafile)
    print(intMtx)
    transMtx  = intMtx.transferMtx
    survivMtx = intMtx.transferAMtx
    sigmaMtx  = intMtx.sigmaMtx
    ydsdyMtx  = intMtx.inelasticityMtx

def getIntMtx(javafile):
    with open(javafile,"rb") as f:
        sobj = f.read()

    intMtx = javaobj.loads(sobj)
    return intMtx

if __name__ == "__main__":
    deserialize()
