# ${renderer.title}

${'##'} Required keywords
%for kw in renderer.required_keywords:
- [${renderer.get_kw_repr(kw)}](#${kw.name})
%endfor

${'##'} Optional keywords
%for kw in renderer.optional_keywords:
- [${renderer.get_kw_repr(kw)}](#${kw.name})
%endfor

${'##'} Keyword documentation

%for groups in (renderer.structured_groups, renderer.unstructured_groups):

%if loop.index==1 and len(groups):
${'###'} Misc. parameters
%endif

%for group in groups:

%if loop.index>0:
##----
<div style='height:20px;'></div>
%endif

<div class='keyword-group' markdown="1">

<%
doc = renderer.get_doc(group)
%>\
${doc}

%if len(doc):
----
%endif

%for kw in group.keywords:

## If this is the second keyword in a group of keywords, and the previous
## keyword had some documentation, then add a separator
%if loop.index>0:
%if len(doc):
----
%endif
%endif

<%
doc = renderer.get_doc(kw).strip()
%>\


## Add some extra margin at the end of the keyword if this is the last keyword
## in the group and there is no documentation (because the margin comes automatically
## in the <p></p> block from the documentation).
%if len(doc) or not loop.last:
${'####'} ${renderer.get_kw_repr(kw)} {:#${kw.name} style='margin:0;'}
%else:
${'####'} ${renderer.get_kw_repr(kw)} {:#${kw.name} style='margin:0; margin-bottom:16px;'}
%endif

%if len(doc):

${doc}

%endif

%endfor ##for for kw in group.keywords

</div>

%endfor ##for group in groups
%endfor ##for groups in (renderer.structured_groups, renderer.unstructured_groups)
