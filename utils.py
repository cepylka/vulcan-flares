import pandas
import numpy


def readSfluxFromFile(filePath, fileIsPlainText, firstRowOnly=True):
    sfluxData = None

    if fileIsPlainText:
        sfluxData = numpy.genfromtxt(
            filePath,
            dtype=float,
            skip_header=1,
            names=["lambda", "flux"]
        )
    else:
        sfluxDataRows = numpy.load(
            filePath,
            allow_pickle=True,
            mmap_mode="r"
        )
        # the pickle might contain several rows, and we only need the first one
        if firstRowOnly:
            fluxAtFirstRow = [
                (lmbd, flx) for lmbd, flx in (
                    sfluxDataRows.iloc[0].items()
                )
            ]
            fluxAtFirstRow_dtype = numpy.dtype("float64, float64")
            sfluxData = numpy.array(fluxAtFirstRow, dtype=fluxAtFirstRow_dtype)
            sfluxData.dtype.names = ["lambda", "flux"]
        else:
            raise NotImplementedError(
                " ".join((
                    "Reading all the rows isn't implemented",
                    "(as it wasn't needed before)"
                ))
            )

    return sfluxData
