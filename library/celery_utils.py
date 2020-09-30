"""
Created on 23/07/2014

From: http://stackoverflow.com/a/23908345

"""
from celery.result import AsyncResult


def store(node):
    id_chain = []
    while node.parent:
        id_chain.append(node.id)
        node = node.parent
    id_chain.append(node.id)
    return id_chain


def restore(id_chain):
    id_chain.reverse()
    last_result = None
    for tid in id_chain:
        result = AsyncResult(tid)
        result.parent = last_result
        last_result = result
    return last_result


def execute_task(task, run_async=False):
    """ execute celery task async or not
        Returns result if not async
    """

    if run_async:
        result = task.apply_async()
    else:
        result = task.apply()
        if not result.successful():
            raise Exception(result.result)

    return result
