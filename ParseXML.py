import re
import sys
fileName = sys.argv[1]

with open(fileName) as infl:
    lines = infl.readlines()

all = ''.join(lines)

chunks = all.split(r'<spectrum_query')
del chunks[0]

pattern = r'([a-zA-Z0-9_]+=)(\"(.*?)\")'
with open('tmp.txt', 'wb') as outfl:
    for chunk in chunks:
        out = []
        chunk = chunk.strip().split('\n')
        chunk = map(lambda x: x.strip(), chunk)
        chunk = filter(lambda x: '=' in x, chunk)
        results = map(lambda x: re.findall(pattern, x), chunk)
        for res in results:
            for r in res:
                out.append(''.join(list(r[0:-1])))
        out = map(lambda x: x.replace('"', ''), out)
        sys.stdout.write('@'.join(out)) #<--- field separator
        sys.stdout.write('\$')           #<--- field chunk
        outfl.writelines('\n'.join(out))
        outfl.write('\n=============================================================\n')
