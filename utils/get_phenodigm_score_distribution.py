import gzip
import optparse
import json
import logging
import csv
import numpy as np
import pandas as pd



LIMIT = 4000000
THRESHOLD = 0.60

def main():
    parser = optparse.OptionParser()
    parser.add_option('-i', '--input', type='string', default=None, dest='filename')
    parser.add_option('-o', '--output', type='string', default=None, dest='outputFilename')
    parser.add_option('-d', '--dataframe', type='string', default=None, dest='dataframeFilename')
    parser.add_option('-s', '--score', type='string', default=None, dest='score_threshold')
    parser.add_option('-f', '--filtered', type='string', default=None, dest='filtered_evidence')


    options, args = parser.parse_args()
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

    logging.info("Read Phenodigm evidence")

    score_threshold = THRESHOLD
    if options.score_threshold:
        score_threshold = float(options.score_threshold)

    scores = np.empty([1, 1], dtype=float)
    score_list = []
    target_list = []
    model_list = []
    disease_list = []
    count = 0
    count_above = 0

    diseases = set()
    targets = set()
    disease_models = set()

    filtered_evidence_fh = open(options.filtered_evidence, mode="wb")

    with gzip.GzipFile(filename=options.filename, mode="rb") as fh:
        for line in fh:
            count+=1
            obj = json.loads(line)
            score = float(obj['evidence']['disease_model_association']['resource_score']['value'])
            if score >= score_threshold:
                #scores = np.append(scores, score)
                score_list.append(score)
                target_list.append(obj['target']['id'])
                model_list.append(obj['evidence']['biological_model']['model_id'])
                disease_list.append(obj['disease']['id'])
                count_above +=1
                filtered_evidence_fh.write(line)
            diseases.add(obj['disease']['id'])
            targets.add(obj['target']['id'])
            disease_models.add(obj['evidence']['biological_model']['model_id'])
            #scores.append(score)
            if count % 10000 == 0:
                print("%i %i"%(count, count_above))
            #if count > LIMIT:
            #    break
    #np.savetxt(options.outputFilename, scores, fmt='%.5f')
    # copy score_list to numpy score

    filtered_evidence_fh.close()

    scores = np.asarray(score_list)

    d = dict(
        score=score_list,
        target=target_list,
        model=model_list,
        disease=disease_list
    )
    df = pd.DataFrame(d)
    df.to_pickle(options.dataframeFilename)
    # save

    print("#scores in np.asarray: %i"%(scores.size))
    print(count_above)
    print("#diseases: %i"%len(diseases))
    print("#targets: %i" % len(targets))
    print("#disease models: %i" % len(disease_models))
    np.save(options.outputFilename, scores)

    # Write results to file
    '''
    with open(options.outputFilename, mode='wb') as zf:
        writer = csv.writer(zf, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['score'])
        for score in scores:
            writer.writerow([score])
    '''

if __name__ == "__main__":
    main()

