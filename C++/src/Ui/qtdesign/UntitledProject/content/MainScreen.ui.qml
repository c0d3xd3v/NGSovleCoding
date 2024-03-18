

/*
This is a UI file (.ui.qml) that is intended to be edited in Qt Design Studio only.
It is supposed to be strictly declarative and only uses a subset of QML. If you edit
this file manually, you might introduce QML code that is not supported by Qt Design Studio.
Check out https://doc.qt.io/qtcreator/creator-quick-ui-forms.html for details on .ui.qml files.
*/
import QtQuick 6.6
import QtQuick.Controls 6.6
import QtQuick.Layouts 6.6

import UntitledProject

//import VTK 9.3
Pane {
    implicitWidth: 600
    implicitHeight: 400

    anchors.fill: parent
    rightPadding: 0
    bottomPadding: 0
    topPadding: 0
    leftPadding: 0
    padding: 0

    SplitView {
        anchors.fill: parent
        anchors.rightMargin: 0
        orientation: Qt.Horizontal

        Pane {
            id: pane
            bottomPadding: 0
            topPadding: 0
            rightPadding: 0
            leftPadding: 0
            padding: 0
            SplitView.preferredWidth: 250
            SplitView.minimumWidth: 250
            SplitView.fillHeight: True

            SwipeView {
                id: swipeView
                anchors.fill: parent

                Mainpage {}

                Item {
                    id: item1
                    width: 200
                    height: 200

                    Loader {
                        id: loader
                        anchors.fill: parent
                        source: ""
                        onSourceChanged: swipeView.currentIndex = 1
                    }
                }
            }
        }

        Pane {
            id: pane1
            SplitView.preferredWidth: 400
            //SplitView.minimumWidth: 200
            SplitView.fillHeight: True


            /*
            VTKRenderWindow {
                id: vtkwindow
                anchors.fill: parent
                // add one or more vtk render items
                property var renderSceneItem: renderItem
                VTKRenderItem {
                    id: renderItem
                    objectName: "ConeView"
                    anchors.fill: vtkwindowr
                    // Provide the handle to the render window
                    renderWindow: vtkwindow
                }
            }
            */
        }
    }
}
