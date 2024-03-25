import QtQuick.Window
import QtQuick.Controls
import QtQuick.Layouts
import QmlVtk 1.0

Window {
    width: 600
    height: 600
    visible: true

    SplitView {
        id: splitView
        anchors.fill: parent

        handle: Rectangle {
                id: handleDelegate
                implicitWidth: 1
                implicitHeight: splitView.height
                color: SplitHandle.pressed ? "#AAAAA"
                    : (SplitHandle.hovered ? Qt.lighter("#AAAAA", 1.5) : "#AAAAA")

                containmentMask: Item {
                    x: (handleDelegate.width - width) / 2
                    width: 10
                    height: splitView.height
                }
            }

        Test {
            //SplitView.fillHeight: true
            //SplitView.fillWidth: true
            SplitView.minimumWidth: 240
        }

        VTKItem {
            objectName: "ConeView"
            //SplitView.fillHeight: true
            //SplitView.fillWidth: true
        }
    }
}
