from pathlib import Path
from fractions import Fraction
import argparse,json
p=argparse.ArgumentParser();p.add_argument("source",type=Path);p.add_argument("output",type=Path);a=p.parse_args()
if a.output.exists():raise FileExistsError("Refuse to overwrite an existing certificate")
c=json.loads(a.source.read_text());out=[];bound=Fraction(0);removed=0
for block in c["core"]["blocks"]:
 rows=[]
 for row in block["factor"]:
  s=sum(map(abs,row))
  if s<1000:removed+=1;bound+=Fraction(s*s,c["core"]["denominator"]**2)
  else:rows.append(row)
 if rows:out.append({**block,"factor":rows})
c["core"]["blocks"]=out
a.output.write_text(json.dumps(c,separators=(",",":")))
print(json.dumps({"removed_rows":removed,"discarded_square_norm_bound_Ha":str(bound),"status":"must independently replay the modified certificate"}))
