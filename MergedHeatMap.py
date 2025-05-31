from utils import *



class MergedHeatMap(object):
    """
    Merge cancer heatmap and lymp heatmap together.
    Input image arrays should be opencv type(BGR).
    """
    def __init__(self, cancerImg, lympImg):
        self.cancerImg = cancerImg
        self.lympImg = lympImg
        self.mergedHeatMap = self.merge()

    def merge(self):
        lympImg = self.lympImg
        cancerArray = self.cancerImg

        up = int(math.ceil(lympImg.shape[0]/cancerArray.shape[0]))

        if up > 1:
            iml_u = np.zeros((cancerArray.shape[0] * up, cancerArray.shape[1] * up), dtype=np.float32)
            for x in range(cancerArray.shape[1]):
                for y in range(cancerArray.shape[0]):
                    iml_u[y * up:(y + 1) * up, x * up:(x + 1) * up] = cancerArray[y, x]

            cancerArray = iml_u.astype(np.float64)

        smooth5 = cancerArray

        smooth5 = cv2.resize(smooth5, (lympImg.shape[1], lympImg.shape[0]), interpolation=cv2.INTER_LINEAR)
        smooth5 = cv2.GaussianBlur(smooth5, (5, 5), 0)

        out = np.zeros((lympImg.shape[0], lympImg.shape[1], 3), dtype=np.uint8)

        for i in range(lympImg.shape[0]):
            for j in range(lympImg.shape[1]):
                out[i, j] = np.array([192,192,192])
                is_tumor = smooth5[i, j] > 0.5
                is_lym = lympImg[i, j] > 0.5
                is_tisue = (smooth5[i, j]!=0) | (lympImg[i, j]!=0)
                # Tissue, Tumor, Lym
                if (not is_tumor) and (not is_lym): # BGR
                    if not is_tisue:
                        out[i, j] = np.array([255,255,255]) #White
                    else:
                        out[i, j] = np.array([192,192,192]) # Grey
                elif is_tumor and (not is_lym):
                    out[i, j] = np.array([0,255,255]) #Yellow
                elif (not is_tumor) and is_lym:
                    out[i, j] = np.array([0,0,200]) #Redish
                else:
                    out[i, j] = np.array([0,0,255]) #Red
        return out
