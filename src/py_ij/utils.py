import imagej

ij = imagej.init()
# Load an image.
image_url = 'https://samples.fiji.sc/new-lenna.jpg'
jimage = ij.io().open(image_url)
# Convert the image from ImageJ to xarray, a package that adds
# labeled datasets to numpy (http://xarray.pydata.org/en/stable/).
image = ij.py.from_java(jimage)

# Display the image (backed by matplotlib).
ij.py.show(image, cmap='gray')