"""Pull the complete list back to two explicit higher-genus covers."""
import json
from points import DEFAULT,classify
from prior.cover_transfer import plane_cover,power_cover


def classify_covers():
    # The inherited algorithms extract exact affine fibres and normalize infinity.
    base=classify(DEFAULT)
    curves=[
        power_cover(DEFAULT,1,2,base),
        plane_cover(DEFAULT,2,'1/2','1/2',base),
    ]
    if [(r['genus'],r['point_count'],r['status']) for r in curves]!=[(5,24,'complete'),(55,4,'complete')]:
        raise ArithmeticError('Transferred cover classification failed')
    return curves


if __name__=='__main__':
    print(json.dumps(classify_covers(),indent=2))
