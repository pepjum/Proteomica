#print 'Dame el file payo!'
import os
import sys

# path = os.getcwd()
path = str(sys.argv[1])
files = os.listdir(path)
files = filter(lambda x: os.path.isfile(os.path.join(path, x)), files)
sqtFiles = filter(lambda x: (x.endswith('sqt')) and ('Corr' not in x), files)


for infile in sqtFiles:
    infl = open(os.path.join(path,infile), 'rb')
    outName = infile.split(r'.')[0] + "_Corr.sqt"
    print("Writing file...{}".format(outName))
    outfl = open(os.path.join(path,outName), 'wb')
    NewValue = '0.0'

    for ln in infl:
        if ln.startswith('S'):
            newLine = ln.split('\t')
            if (newLine[7] == '') and (newLine[8] == ''):
                newLine[7] = NewValue
                newLine[8] = NewValue
            outfl.write('\t'.join(newLine))
        else:
            outfl.write(ln)

    outfl.close()
    infl.close()

print('Done!!')
