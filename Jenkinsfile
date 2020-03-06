pipeline {
    agent any

    stages {
        stage('Checkout') { // Checkout (git clone ...) the projects repository
            steps {
            checkout scm
            }
        }
        stage('Build') {
            steps {
                echo 'Building..'
                sh '/home/tyler/anaconda3/bin/python -m venv env'
                sh 'source env/bin/activate'
                sh 'python -m pip install -r requirements.txt'
            }
        }
        stage('Test') {
            steps {
                echo 'Testing..'
                sh '''
                python vase.py -h
                '''
            }
        }
    }
}
