import io, json, tempfile, unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch
import lambda_pool

class LambdaPoolTests(unittest.TestCase):
    def test_launch_body_and_redaction(self):
        with tempfile.TemporaryDirectory() as d:
            old = lambda_pool.REGISTRY; lambda_pool.REGISTRY = Path(d)/"pool.json"
            with patch.object(lambda_pool, "_json", side_effect=[{"data":{"gpu_1x_a10":{"instance_type":{"price_cents_per_hour":129},"regions_with_capacity_available":[{"name":"us-west-1"}]}}},{"data":{"instance_ids":["i-1"]}}]) as call:
                lambda_pool.launch(type("A", (), {"instance_type":"A10","name":"test","region":"us-west-1"})())
                self.assertEqual(call.call_args.args[0:2], ("POST", "/instance-operations/launch")); self.assertEqual(call.call_args.args[2]["quantity"], 1); self.assertEqual(call.call_args.args[2]["ssh_key_names"], ["aiden-mac"])
            lambda_pool.REGISTRY = old
    def test_unknown_termination_refused(self):
        with tempfile.TemporaryDirectory() as d:
            old = lambda_pool.REGISTRY; lambda_pool.REGISTRY = Path(d)/"pool.json"
            with self.assertRaises(SystemExit): lambda_pool.terminate(type("A", (), {"id":"unknown"})())
            lambda_pool.REGISTRY = old
    def test_safe_redacts_jupyter(self):
        self.assertNotIn("token", json.dumps(lambda_pool._safe({"url":"https://x/jupyter/token"})))
    def test_status_main_allowlist(self):
        output=io.StringIO()
        with patch.object(lambda_pool,'_json',return_value={'data':[{'id':'i','status':'active','jupyter_url':'secret-url','secret':'secret-value'}]}), redirect_stdout(output):
            lambda_pool.main(['status'])
        self.assertEqual(json.loads(output.getvalue()),[{'id':'i','status':'active'}])
    def test_price_and_region_refuse_without_post(self):
        for price,region in [(130,'us-west-1'),(129,'no-capacity')]:
            data={'data':{'gpu_1x_a10':{'instance_type':{'price_cents_per_hour':price},'regions_with_capacity_available':[{'name':'us-west-1'}]}}}
            with patch.object(lambda_pool,'_json',return_value=data) as call:
                with self.assertRaises(SystemExit):lambda_pool.main(['launch','--name','test','--region',region])
                self.assertEqual(call.call_count,1)
                self.assertEqual(call.call_args.args[0],'GET')
    def test_register_actual_shape_idempotent(self):
        with tempfile.TemporaryDirectory() as d, patch.object(lambda_pool,'REGISTRY',Path(d)/'pool.json'):
            data={'data':[{'id':'i','instance_type':{'name':'gpu_1x_a10'},'region':{'name':'us-west-1'}}]}
            with patch.object(lambda_pool,'_json',return_value=data):
                lambda_pool.main(['register','i']);lambda_pool.main(['register','i'])
            self.assertEqual(len(lambda_pool._read_registry()['instances']),1)

if __name__ == "__main__": unittest.main()
