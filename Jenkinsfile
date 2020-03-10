pipeline {
    agent any

    stages {
        stage('Checkout') {
            steps {
            checkout scm
            }
        }
        stage('Build Environment') {
            steps {
                sh '''
                python -m venv env
                source env/bin/activate
                python -m pip install -r requirements.txt
                python -m pip install pyenchant
                python -m pip install pylint
                '''
            }
        }

        stage('Pylint') {
            steps {
                sh '''
                source env/bin/activate
                python -m pylint --exit-zero -f parseable -d W1202 --ignore=config,deprecated,docs,env,tests,VaSeEval,VaSeUtils,vaseutils.py ./*py > pylint.report
                '''
            }
            post {
                always{
                    recordIssues(
                        tools: [pyLint(pattern: 'pylint.report')],
                        unstableTotalAll: 100,
                        failedTotalAll: 120,
                        unstableTotalNormal: 10,
                        failedTotalNormal: 15,
                        failedTotalHigh: 1
                    )
                }
            }
        }
        stage('Help Menues') {
            steps {
                sh '''
                source env/bin/activate
                python vase.py -h
                python vase.py BuildSpikeIns -h
                python vase.py AssembleValidationSet -h
                python vase.py BuildValidationSet -h
                '''
            }
        }
    }
}
