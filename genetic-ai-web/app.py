
"""
@author__ = "Juan Francisco Illan"
@license__ = "GPL"
@version__ = "1.0.1"
@email__ = "juanfrancisco.illan@gmail.com"
"""

from flask import Flask, request, render_template
from flask_debugtoolbar import DebugToolbarExtension

app = Flask(__name__)
app.config['SECRET_KEY'] = "lailolailo"
app.debug = True
toolbar = DebugToolbarExtension(app)

# Controller to index
@app.route('/')
def home():
	return render_template('index.html')
	

# Controller to blast access
@app.route('/about', methods=["GET"])
def about():

	return render_template('about.html')

# startup HTTP web service
if __name__ == '__main__':
	print("Running service ...")
	app.run(host='0.0.0.0')

