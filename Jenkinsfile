pipeline {
    agent any

    stages {
        stage('Build') {
            steps {
                echo 'Building..'
                /home/tyler/anaconda3/bin/python -m venv env -m venv env
                source env/bin/activate
                python -m pip install -r requirements.txt
            }
        }
        stage('Test') {
            steps {
                echo 'Testing..'
                python vase.py -h
            }
        }
    }
}
