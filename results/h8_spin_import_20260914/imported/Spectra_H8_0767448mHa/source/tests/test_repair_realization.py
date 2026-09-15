import unittest,copy
from fractions import Fraction as F
from experiments.repair_realization import learn,verify,verify_against_records,state_after,mv,step


class RealizationTests(unittest.TestCase):
    def evidence(self):
        # An independent 3-state shift system: delayed, then decaying output.
        z=[F(0)]*3;records=[];history=[]
        for t in range(8):
            u=F(t==0);history.append(u)
            z=[u,z[0],z[1]+z[2]/2]
            records.append(dict(controls=list(history),output_intervals=[[z[2],z[2]]],reset_confirmed=True))
        return dict(prior='zero_reset_strictly_proper_LTI',maximum_dimension=4,output_ports=1,records=records)

    def test_discovers_delayed_rank_and_updates(self):
        data=self.evidence();m=learn(data)['model'];self.assertEqual(m['rank'],3)
        self.assertEqual(verify(m)['status'],'exact_finite_order_certificate')
        z=[F(0)]*3;u=[F(1,3),F(1,2),F(0),F(1,5),F(0)]
        for v in u:z=[v,z[0],z[1]+z[2]/2]
        self.assertEqual(mv(m['output'],state_after(m,u,reset_confirmed=True)),[z[2]])

    def test_insufficient_data_requests_experiment(self):
        data=self.evidence();data['records']=data['records'][:2]
        answer=learn(data);self.assertEqual(answer['status'],'needs_evidence')
        self.assertEqual(set(answer['identified_markov_vectors']),{'0','1'})
        self.assertEqual(len(answer['experiment']['controls']),8)

    def test_uncertainty_cannot_be_rounded_into_proof(self):
        data=self.evidence();data['records'][0]['output_intervals'][0][1]=F(1,1000000)
        with self.assertRaisesRegex(ValueError,'uncertain'):learn(data)

    def test_tamper_initialization_and_interfaces(self):
        m=learn(self.evidence())['model'];bad=copy.deepcopy(m);bad['input'][0]+=1
        with self.assertRaises(ValueError):verify(bad)
        with self.assertRaisesRegex(ValueError,'recorded history'):state_after(m,None)
        with self.assertRaises(ValueError):step(m,m['reset'],F(2))
        data=self.evidence();data['records'][0]['reset_confirmed']=False
        with self.assertRaisesRegex(ValueError,'initialization'):learn(data)

    def test_missing_order_bound_and_inconsistent_records(self):
        data=self.evidence();data['maximum_dimension']=0
        with self.assertRaises(ValueError):learn(data)
        data=self.evidence();r=copy.deepcopy(data['records'][0]);r['output_intervals']=[[F(1),F(1)]]
        data['records'].append(r)
        with self.assertRaisesRegex(ValueError,'inconsistent'):learn(data)

    def test_zero_system(self):
        data=self.evidence()
        for r in data['records']:r['output_intervals']=[[0,0]]
        m=learn(data)['model'];self.assertEqual(m['rank'],0)
        self.assertEqual(mv(m['output'],state_after(m,[1,0],reset_confirmed=True)),[0])

    def test_self_consistent_fabrication_cannot_replace_evidence(self):
        data=self.evidence();m=learn(data)['model']
        self.assertEqual(verify_against_records(m,data)['status'],'verified_against_observable_records')
        bad=copy.deepcopy(m);bad['input']=[2*x for x in bad['input']]
        bad['observed_moments']=[[2*x for x in row] for row in bad['observed_moments']]
        verify(bad)  # Internally consistent is insufficient.
        with self.assertRaisesRegex(ValueError,'observable evidence'):verify_against_records(bad,data)

    def test_malformed_evidence_has_explicit_refusal(self):
        for data in ({}, [], None):
            with self.subTest(data=data):
                with self.assertRaisesRegex(ValueError,'evidence mapping'):learn(data)
        for record in ({}, [], None):
            data=self.evidence();data['records'][0]=record
            with self.subTest(record=record):
                with self.assertRaisesRegex(ValueError,'record mapping'):learn(data)


if __name__=='__main__':unittest.main()
