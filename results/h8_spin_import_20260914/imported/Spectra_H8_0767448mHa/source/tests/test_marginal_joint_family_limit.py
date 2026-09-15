import unittest
from fractions import Fraction as F
from experiments.marginal_joint_family_limit import replay


class JointFamilyLimitTests(unittest.TestCase):
    def certificate(self):
        return {'kind':'joint_projector_family_limit_v1','half_vector':{0x999:1,0x666:1},
                'charged_vector':{62:1,3008:1},'ratio':'1/2','theta_half':'1/2',
                'theta_joint':'3/4','mixture':[{'weight':'1','vector':{0:1}}]}

    def test_vacuum_and_full_mixture_independent_expectation(self):
        # Vacuum and fully occupied determinants both have centered onsite
        # energy10, density energy5/2, and zero hopping expectation for EVERY
        # mean-correct profile, giving the deliberately loose cap5/2.
        c=self.certificate();c['mixture']=[{'weight':'1/3','vector':{0:7}},
                                         {'weight':'2/3','vector':{4095:11}}]
        r=replay(c)
        self.assertEqual(F(r['periodic_family_upper']),F(5,2))
        self.assertEqual(r['profile_gradients'],['0']*6)
        self.assertEqual(r['half_expectation'],'0')
        self.assertEqual(r['joint_expectation'],'0')
        self.assertEqual(r['mixture_trace'],'1')

    def test_refuses_trace_positivity_profile_and_fidelity_errors(self):
        for mixture in ([{'weight':'1/2','vector':{0:1}}],
                        [{'weight':'-1','vector':{0:1}},{'weight':'2','vector':{4095:1}}],
                        [{'weight':'1','vector':{0:0}}],
                        [{'weight':'1','vector':{1:1}}],
                        [{'weight':'1','vector':{0x999:1,0x666:1}}]):
            c=self.certificate();c['mixture']=mixture
            with self.subTest(mixture=mixture),self.assertRaises(ValueError):replay(c)
        c=self.certificate();c['theta_joint']='0';c['mixture']=[{'weight':'1','vector':{62:1,3008:1}}]
        with self.assertRaisesRegex(ValueError,'fidelity'):replay(c)
        c=self.certificate();c['proposed_periodic_family_upper']='0'
        with self.assertRaisesRegex(ValueError,'disagrees'):replay(c)

    def test_refuses_malformed_sources_and_claims(self):
        for field,value in [('ratio','0'),('theta_half','2'),('theta_joint','2'),
                            ('half_vector',{0x555:1}),('charged_vector',{62:1}),('mixture',[]),
                            ('kind','ground_energy_upper')]:
            c=self.certificate();c[field]=value
            with self.subTest(field=field),self.assertRaises(ValueError):replay(c)


if __name__=='__main__':unittest.main()
