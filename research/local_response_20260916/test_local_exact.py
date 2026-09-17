import copy
from fractions import Fraction as F
import unittest
from research.local_response_20260916 import local_exact as c


def payload():
    return {'kind':'local_hubbard_boundary_v1','rungs':4,'U':'8','t':'1',
            'messages':[c.zero(16),c.zero(16)],'lower_shifts':['-10']*3}


def mixed_dual():
    p=payload();p['dual_marginals']=[]
    for _ in range(3):
        blocks={}
        for key,labels in c.sector_labels(4).items():
            blocks[f'{key[0]},{key[1]}']=[
                [F(int(i==j),256) for j in labels] for i in labels]
        p['dual_marginals'].append(blocks)
    return p


class LocalExactTests(unittest.TestCase):
    def test_rung_and_all_leg_signs(self):
        h=c.rung_matrix()
        self.assertEqual(h[1][4],-1)
        self.assertEqual(h[9][3],-1)
        self.assertEqual(h[6][3],1)
        for blocks in c.local_blocks([c.zero(16),c.zero(16)]):
            for src,dst in ((0,4),(1,5),(2,6),(3,7)):
                a,b=1<<src,1<<dst
                labels,m=blocks[c.spin_counts(a,4)]
                self.assertEqual(m[labels.index(a)][labels.index(b)],-1)

    def test_local_coverage_and_acceptance(self):
        labels=c.sector_labels(4)
        self.assertEqual(sum(map(len,labels.values())),256)
        self.assertEqual(max(map(len,labels.values())),36)
        rec=c.verify(payload())
        self.assertEqual(rec['energy_lower_over_t'],'-30')
        self.assertEqual(rec['global_determinants_enumerated'],0)

    def test_valid_offdiagonal_and_refusals(self):
        p=payload();p['messages'][0]=c.rung_matrix()
        c.verify(p)
        for mode in ('asymmetry','wrong_charge','bad_lower','wrong_model'):
            bad=payload()
            if mode=='asymmetry':bad['messages'][0][1][4]=1
            if mode=='wrong_charge':bad['messages'][0][0][1]=bad['messages'][0][1][0]=1
            if mode=='bad_lower':bad['lower_shifts']=['100']*3
            if mode=='wrong_model':bad['t']='2'
            with self.assertRaises(ValueError,msg=mode):c.verify(bad)

    def test_dual_and_refusals(self):
        good=mixed_dual()
        rec=c.verify_dual(good)
        self.assertEqual(F(rec['family_lower_ceiling_over_t']),16)
        for mode in ('trace','negative','mismatch'):
            bad=copy.deepcopy(good)
            a=bad['dual_marginals'][0]
            a['0,0'][0][0]+=F(1,512)
            if mode in ('negative','mismatch'):
                a['4,4'][0][0]-=F(1,512) if mode=='mismatch' else 1
            with self.assertRaises(ValueError,msg=mode):c.verify_dual(bad)

    def test_complete_identity_independently_on_4900_state_sector(self):
        # Charged test ONLY. Accepting local checker never calls this model.
        from research.constructive_response_20260916.interacting_probe import rational_rows
        labels,expected=rational_rows(4,F(1));index={x:i for i,x in enumerate(labels)}
        messages=[c.rung_matrix(),[[2*x for x in row] for row in c.rung_matrix()]]
        blocks=c.local_blocks(messages)
        actions=[]
        for block in blocks:
            table={}
            for local_labels,m in block.values():
                for j,state in enumerate(local_labels):
                    table[state]={target:m[i][j] for i,target in enumerate(local_labels) if m[i][j]}
            actions.append(table)
        for col,state in enumerate(labels):
            row={}
            for bi,table in enumerate(actions):
                shift=4*bi;mask=255<<shift;local=(state>>shift)&255
                for target,value in table[local].items():
                    global_target=(state&~mask)|(target<<shift)
                    i=index[global_target];row[i]=row.get(i,F(0))+value
            row={i:v for i,v in row.items() if v}
            # Both models are symmetric; expected[col] is a row map.
            self.assertEqual(row,expected[col])

if __name__=='__main__':unittest.main()
