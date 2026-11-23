# The original Part XXX solver

The atlas can run the archived Part XXX decision procedure directly on a supplied quartic. The ten Python modules in `vendor/crest` are byte-for-byte copies of the modules used by that procedure in `crest_root_calculus_integrated_part_xxx.zip`. [Source hashes](source-hashes.json) identify the archive and each module. The adapter converts a homogeneous quartic to a list of integer monomials, checks smoothness, and tries the supplied primes in order.

Install the optional scientific dependencies from this folder with `python -m pip install -r requirements.txt`. Then run:

```sh
python ../app.py ../original-slice.json
python ../app.py ../original-C-star.json
python ../app.py ../original-nonlift.json
python replay.py
```

The blind twelve-row sample comes from the same equation list as the portable atlas. The original engine certifies zero geometric degree-three classes for eleven rows. The last row, **23640**, requires a larger scan at the eligible primes. It stays open under the shipped field-size budget. The bounds meet at **one** for the C-star control (prime 5) and at **zero** for the non-lifting control (prime 13). The positive receipt contains its Chow polynomial, field and flags. See the [catalogue](catalogue.json), [positive result](C-star.json) and [non-lift result](nonlift.json).

To test your own curve, replace the `id` and `quartic` values in a Part XXX case. The local browser and command line share the same adapter. `max_field_size` gates expensive extension-field scans and `max_precision` limits lifting. An unresolved row carries the reason. Inspection of the exact algorithm itself can start at [decide.py](vendor/crest/scan_completeness_xxx/decision/decide.py). The archive’s standalone driver accepts the original monomial-list format and offers its original unrestricted geometric scope.

Run `python audit_powers.py` here for a separate characteristic-zero audit of the portable solver’s seven Frobenius exclusions. It checks 392 powered sextics independently with SymPy.
