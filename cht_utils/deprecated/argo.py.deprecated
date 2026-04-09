from hera.workflows import Workflow, WorkflowStatus, Task
from hera.workflows.models import WorkflowTemplateRef
from hera.shared import GlobalConfig
import time

class Argo:

    def __init__(self, host, workflow_name, wait=True):
        self.host = host
        self.workflow_name = workflow_name
        self.wait = wait

    def submit_job(self, **kwargs):
        GlobalConfig.namespace = "argo"
        GlobalConfig.host = self.host
        GlobalConfig.verify_ssl = False
        wt_ref = WorkflowTemplateRef(name=self.workflow_name,
                                     cluster_scope=False)
        w = Workflow(
            generate_name=self.workflow_name+"-",
            workflow_template_ref=wt_ref,
            arguments=kwargs
        )
        w.create()
        if self.wait:
            time.sleep(1) 
            w.wait()
        return w
