# The IDS *path* syntax to identify a part of an IDS
​
## Scope
​
This document describes the syntax of the ***path*** that identifies a *part* of an IDS.
In this context a part can be defined by a given node (or node element in the case of an array of structure) 
with all its leaves, or by a single leaf of the IDS structure. As a reminder, nodes are structures or arrays 
of structure and leaves are quantities. 
​
Parts that are represented by arrays (either arrays of structure or quantities) can be identified
as a whole (all the elements of the array), as a range (contiguous elements within two bounds) or 
as a single element.

This document does not cover the identification of the IDS itself, please refer to fragment section 
of the It can be used standalone or in complement of the [IMAS URI syntax](IMAS-URI-scheme.md) document
for that purpose.
​
​
## Syntax
​
An IDS *path* is represented by a string that respects the following rules:
​
- a slash **`/`** separates nested structures, the leading one indicates the root of an IDS
- scalar quantities as well as structures are given by their name
- arrays (quantities and arrays of structure) have their names followed by round brackets **`()`**, this is mandatory unless this array is in last position in the path 
- when brackets are added, a **range specifier** is mandatory inside, e.g in 1D:
	+ `(:)` indicates all elements of the array
	+ `(i)` indicates the i-th element of the array, 1 being the first element (same as Fortran arrays)
	+ `(i:j)` indicates all elements from the i-th to the j-th (both included)
- for nD arrays, a list of n range specifiers separated by commas must be given, e.g in 3D 
	+ `(:,:,1)` represents a 2D slab for the first element in the last dimension (usually the time, unless the quantity is 
	not time dependent or has one time dependent array of structure as ancester)
​
​
## Use cases and examples
​
An IDS *path* can be used in the following cases: 

- alone, it can be used as a value of some fields (string) in IDSs or as an argument of a `partial_get` operation from a data-access library;
- as part of a URI fragment to identify a specific subset of a given IDS in a given data-entry.
​
### *idspath* or *fragment* as `path` field in IDSs

In several places in the Data-Dictionary, specific string fields of IDS allow to make a reference to other fields of the same IDS. 
This is the case for instance of the field `ids_properties/provenance/node(:)/path`. In this field, an idspath is expected point to 
a subset of the current IDS.

In other cases, the reference `path` field allows to point to different IDS of the same data-entry (e.g in `grid_ggd(:)/path`). 
In such a case, the string following the syntax of idspath can point to a subset of the same IDS, or a subset of another IDS.
In the latter case, `path` can contain a *same-document* URI fragment. Refer to [URI syntax](IMAS-URI-scheme.md) documentation for more information.

### *idspath* in URI fragment

When the *idspath* is given as a part of an IMAS URI, it is always given in complement information on the IDS (and optionnally its occurrence), 
and follows the syntax `ids[:occurrence][/idspath]` described in the [URI syntax](IMAS-URI-scheme.md) documentation. 

### *idspath* as an argument of a `partial_get` operation
​
When given as an argument of the `partial_get` function of the Access-Layer, please refer to the [Access-Layer User Guide](https://user.iter.org/?uid=YSQENW&action=get_document) 
to get more information about possible limitations of the supported languages and allowed comintations of range specifiers.
​
​