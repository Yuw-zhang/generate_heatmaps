from utils import *
from FourPanelImage import FourPanelImage
from HeatMap import HeatMap

# these folders will be replaced by paramaters
svs_fol = '/data10/shared/tcga_all/brca'
cancer_fol = './prediction_brca/tumor'
til_fol = './prediction_brca/til'
output_pred = 'brca_demo_heatmap'

prefix = "prediction-"
wsi_extension = ".svs"
skip_first_line_pred = False
is_cancer_wsiID_same_til_wsiID = True

fns = [fn.split('prediction-')[-1] for fn in os.listdir(til_fol) if fn.startswith('prediction-') and not fn.endswith('low_res')]
til_wsiID_map = collections.defaultdict(str)
for fn in fns:
    til_wsiID_map[fn.split('.')[0]] = fn


def checkFileExisting(wsiId):
    til_wsiID = wsiId
    if not is_cancer_wsiID_same_til_wsiID:
        til_wsiID = til_wsiID_map[wsiId]  # if cancer id is different from til slide id

    allPath = [
        os.path.join(svs_fol, wsiId + wsi_extension), # svsPath
        os.path.join(cancer_fol, prefix + wsiId), # predPath
        os.path.join(til_fol, 'prediction-' + til_wsiID), # tilPath_pred
    ]
    ans = True
    for path in allPath:
        ans = os.path.exists(path)
        if not ans:
            print(path, "does not exit!")
            break
    return ans


def gen1Image(fn):
    wsiId = fn[len(prefix):]
    if not checkFileExisting(wsiId):
        return

    oslide = openslide.OpenSlide(os.path.join(svs_fol, wsiId + wsi_extension))
    til_wsiID = wsiId
    if not is_cancer_wsiID_same_til_wsiID:
        til_wsiID = til_wsiID_map[wsiId]     # if cancer id is different from til slide id

    til_heatmap = HeatMap(til_fol, skip_first_line_pred=False)
    til_heatmap.setWidthHeightByOSlide(oslide)
    til_map = til_heatmap.getHeatMapByID(til_wsiID)

    cancer_heatmap = HeatMap(cancer_fol, skip_first_line_pred=False)
    cancer_heatmap.setWidthHeightByOSlide(oslide)
    cancer_map = cancer_heatmap.getHeatMapByID(wsiId)

    img = FourPanelImage(oslide, cancer_map, til_map,
                       os.path.join(output_pred, wsiId+".png"))
    img.saveImg()
    print(wsiId)


def main(parallel_processing = 0):
    if not os.path.isdir(output_pred):
        os.makedirs(output_pred)
    print('In main, prefix: ', prefix)

    grades_prediction_fns = [f for f in os.listdir(cancer_fol) if f.startswith(prefix) and not f.endswith('low_res')]

    if len(grades_prediction_fns) == 0:
        print("In main: No valid file!")
        return

    if parallel_processing in {0, 1}:
        for i, fn in enumerate(grades_prediction_fns):
            print(i, fn)
            try:
                gen1Image(fn)
            except:
                pass
    else:
        num_of_cores = multiprocessing.cpu_count() - 2
        if parallel_processing > 0:
            num_of_cores = min(num_of_cores, parallel_processing)

        print("Using multiprocessing, num_cores: ", num_of_cores)
        p = multiprocessing.Pool(num_of_cores)
        p.map(gen1Image, grades_prediction_fns)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        parallel_processing = 0
    else:
        parallel_processing = int(sys.argv[1])

    print("***********************************************")
    print("Usage: python main.py 0/1/4/-1")
    print("0/1: not using parallel processing")
    print("any number larger than 1, N, using N cores in parallel processing")
    print("-1: use all available cores in parallel processing, left 2 cores for others")
    print("***********************************************\n")

    main(parallel_processing = parallel_processing)




