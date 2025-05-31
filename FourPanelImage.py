from MergedHeatMap import MergedHeatMap
from utils import  *


class FourPanelImage(object):
    def __init__(self, oslide, cancerMap, tilMap, savePath):
        self.oslide = oslide
        self.cancerMap = cancerMap
        self.tilMap = tilMap
        self.savePath = savePath

    def saveImg(self):
        shape = (self.cancerMap.shape[1] * 2, self.cancerMap.shape[0] * 2)
        thumbnail = self.oslide.get_thumbnail(shape)

        mergedMap = MergedHeatMap(self.cancerMap, self.tilMap).mergedHeatMap
        mergedMap = mergedMap[:, :, [2, 1, 0]]  # convert to rgb
        mergedMap = np.rot90(mergedMap, k=3)
        mergedMap = np.fliplr(mergedMap)

        cancerSmoothImg = cv2.GaussianBlur(self.cancerMap, (5, 5), 0)
        cancerSmoothImg = np.rot90(cancerSmoothImg, k=3)
        cancerSmoothImg = np.fliplr(cancerSmoothImg)

        TilMapSmoothImg = cv2.GaussianBlur(self.tilMap, (5, 5), 0)
        TilMapSmoothImg = np.rot90(TilMapSmoothImg, k=3)
        TilMapSmoothImg = np.fliplr(TilMapSmoothImg)

        aspect = self.cancerMap.shape[1] / self.cancerMap.shape[0]

        width = 6.4 * 1
        mpl.rcParams["figure.figsize"] = [width * 1.10, width * aspect]
        mpl.rcParams["figure.dpi"] = 300

        if aspect > 1:
            hspace = 0.04
            wspace = hspace * aspect
        else:
            wspace = 0.04
            hspace = wspace / aspect + 0.16# / aspect
        # * aspect
        fig2, axarr = plt.subplots(2, 2, gridspec_kw={'wspace': 0.28, 'hspace': hspace})

        caxarr = []
        for r in range(2):
            for c in range(2):
                divider = make_axes_locatable(axarr[r, c])
                cax = divider.append_axes("right", size="5%", pad=0)
                caxarr.append(cax)
        caxarr = np.array(caxarr).reshape(2, 2)

        for x in [0, 1]:
            for y in [0, 1]:
                axarr[x, y].tick_params(labelbottom=False, labelleft=False, bottom=False, left=False)
                caxarr[x, y].tick_params(labelbottom=False, labelleft=False, bottom=False, left=False)
                caxarr[x, y].axis("off")


        axarr[0, 0].imshow(thumbnail)
        axarr[0, 0].set_title('H&E', fontsize=8)

        # heat map of cancer
        masked_cancer = np.where(cancerSmoothImg == 0, np.nan, cancerSmoothImg)
        cancerIm = axarr[0, 1].imshow(masked_cancer, cmap='jet', vmax=1.0, vmin=0.0)

        axins = inset_axes(caxarr[0, 1],
                           width="50%",  # width = 5% of parent_bbox width
                           height="60%",  # height : 50%
                           loc='lower left',
                           bbox_to_anchor=(0.2, 0, 1, 1),
                           bbox_transform=caxarr[0, 1].transAxes,
                           borderpad=0,
                           )

        cb = fig2.colorbar(cancerIm, cax=axins)
        cb.ax.tick_params(labelsize='xx-small')
        caxarr[0, 1].axis("off")
        axarr[0, 1].set_title('Tumor Spatial Probability Map', fontsize=8)


        # heat map of TIL
        masked_til = np.where(TilMapSmoothImg == 0, np.nan, TilMapSmoothImg)
        TilIm = axarr[1, 0].imshow(masked_til, cmap='jet', vmax=1.0, vmin=0.0)

        axins = inset_axes(caxarr[1, 0],
                           width="50%",  # width = 5% of parent_bbox width
                           height="60%",  # height : 50%
                           loc='lower left',
                           bbox_to_anchor=(0.2, 0, 1, 1),
                           bbox_transform=caxarr[1, 0].transAxes,
                           borderpad=0,
                           )
        cb = fig2.colorbar(TilIm, cax=axins)
        cb.ax.tick_params(labelsize='xx-small')
        caxarr[1, 0].axis("off")
        axarr[1, 0].set_title('Lymphocyte Spatial Probability Map', fontsize=8)


        # megered map
        axarr[1, 1].imshow(mergedMap)
        colors_classification = ['red', 'yellow', 'gray']
        labels_classification = ['Lymph VGG16', 'Tumor', 'Background']
        legend_patches_classification = [Rectangle((0, 0), width=0.01, height=0.01, color=icolor, label=label, lw=0)
                                         for icolor, label in zip(colors_classification, labels_classification)]
        caxarr[1, 1].legend(handles=legend_patches_classification,
                            facecolor=None,  # "white",
                            edgecolor=None,
                            fancybox=False,
                            bbox_to_anchor=(-0.1, 0),
                            loc='lower left',
                            fontsize='x-small',
                            shadow=False,
                            framealpha=0.,
                            borderpad=0)
        axarr[1, 1].set_title('Tumor-TILs Map', fontsize=8)
        plt.savefig(self.savePath, bbox_inches='tight')
        plt.close()
