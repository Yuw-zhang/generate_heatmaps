from PredictionFile import PredictionFile
from utils import *


class HeatMap(object):

    def __init__(self, rootFolder, skip_first_line_pred):
        self.rootFolder = rootFolder
        self.width = 0
        self.height = 0
        self.skip_first_line_pred = skip_first_line_pred

    def setWidthHeightByOSlide(self, slide):
        self.width = slide.dimensions[0]
        self.height = slide.dimensions[1]

    def getHeatMapByID(self, slideId, prefix = 'prediction-'):
        """
        Args:
            slideId(str): id of svs file, like 'TCGA-3C-AALI-01Z-00-DX1'.
        """
        predictionFileName = prefix + slideId
        self.heatmap = self.getHeatMap(os.path.join(self.rootFolder, predictionFileName))
        return self.heatmap

    def getHeatMap(self, predPath):
        """
        Args:
            predPath(str): must be full path.
        """
        predictionFile = PredictionFile(predPath, self.skip_first_line_pred)
        predictionFile.setWidthHeight(self.width, self.height)
        pred, necr, patch_size = predictionFile.get_pred_and_necr()
#         print("pred.shape:",pred.shape)

        return pred
