from django.template import Library

from snpdb.models import DataState

register = Library()


@register.inclusion_tag("seqauto/tags/seqauto_messages.html")
def seqauto_messages(messages):
    return {"messages": messages}


@register.inclusion_tag("seqauto/tags/record_data_state_helper.html")
def record_data_state_helper(record):
    """ Takes a SeqAuto record object (ie something that has a job script attached """

    pbs = None
    pbs_out = None
    try:
        pbs = record.jobscript_set.all().order_by("pk").last()
        if pbs and record.data_state == DataState.ERROR:
            with open(pbs.out_file) as f:
                pbs_out = f.read()
    except AttributeError:
        pass

    return {'record': record,
            'pbs': pbs,
            'pbs_out': pbs_out}
